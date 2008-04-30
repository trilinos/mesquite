/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TrustRegion.cpp
 *  \brief Port Todd Munson's trust region solver to Mesquite
 *  \author Jason Kraftcheck (Mesquite Port)
 */

#include "Mesquite.hpp"
#include "TrustRegion.hpp"
#include "MsqDebug.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

#define USE_FN_PC1    // Use 1st preconditioner from Todd's code
                      // (alternate is whatever is in MsqHessian already)
#undef  DO_STEEP_DESC // Jason's apparently broken hack to fall back to
                      // steepest descent search direction

namespace Mesquite {

msq_std::string TrustRegion::get_name() const { return "TrustRegion"; }

PatchSet* TrustRegion::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

TrustRegion::TrustRegion( ObjectiveFunction* of, bool Nash )
  : VertexMover( of, Nash ), PatchSetUser(true), mMemento(0)
{}

TrustRegion::~TrustRegion()
{
  delete mMemento;
  mMemento = 0;
}


void TrustRegion::initialize( PatchData& pd, MsqError& err )
{
  mMemento = pd.create_vertices_memento( err );
  MSQ_CHKERR(err);
}

void TrustRegion::initialize_mesh_iteration( PatchData& /*pd*/, MsqError& /*err*/ )
{ }

void TrustRegion::terminate_mesh_iteration( PatchData& /*pd*/, MsqError& /*err*/ )
{ }

void TrustRegion::cleanup()
{ 
  delete mMemento;
  mMemento = 0;
}

static inline double norm( const Vector3D* vect, size_t nn )
  { return length_squared( vect, nn ); }

static inline void negate( Vector3D* out, const Vector3D* in, size_t nn )
{
  for (size_t i = 0; i < nn; ++i)
    out[i] = -in[i];
}

static inline void axpy(Vector3D *r, const Vector3D *x, double c, 
                        const Vector3D *y, size_t nn)
{
  for (size_t i = 0; i < nn; ++i) 
    r[i] = x[i] + c*y[i];
}

static inline void maxpy(Vector3D *r, const Vector3D *x, double c, 
                         const Vector3D *y, size_t nn)
{
  for (size_t i = 0; i < nn; ++i)
    r[i] = c*y[i] - x[i];
}

static inline void matmul( Vector3D* result, const MsqHessian& hess, const Vector3D* x )
{
  MsqError err;
  axpy( result, hess.size(), hess, x, hess.size(), 0, 0, err );
}

void TrustRegion::compute_preconditioner( MsqError& err )
{
#ifndef USE_FN_PC1
  mHessian.calculate_preconditioner(err);
#else
  preCond.resize( mHess.size() );
  for (size_t i = 0; i < mHess.size(); ++i) {
    const Matrix3D& m = *mHess.get_block(i,i);
    preCond[i] = 1.0 / (m[0][0] + m[1][1] + m[2][2]);
    if (preCond[i] < 0)
      preCond[i] = 1;
  }
#endif
}

void TrustRegion::apply_preconditioner( Vector3D* z, Vector3D* r, MsqError& err )
{
#ifndef USE_FN_PC1
  mHessian.apply_preconditioner( z, r, err );
#else
  for (size_t i = 0; i < preCond.size(); ++i) 
    z[i] = preCond[i] * r[i];
#endif
}

void TrustRegion::optimize_vertex_positions( PatchData& pd, MsqError& err )
{
  TerminationCriterion& term = *get_inner_termination_criterion();
  OFEvaluator& func = get_objective_function_evaluator();
  
  const double cg_tol = 1e-2;
  const double eta_1  = 0.01;
  const double eta_2  = 0.90;
  const double tr_incr = 10;
  const double tr_decr_def = 0.25;
  const double tr_decr_undef = 0.25;
  const double tr_num_tol = 1e-6;
  const int max_cg_iter = 10000;
  
  double radius = 1000;		/* delta*delta */

  Vector3D *w;
  Vector3D *z;
  Vector3D *d;
  Vector3D *p;
  Vector3D *r;

  double norm_r, norm_g;
  double alpha, beta, kappa;
  double rz, rzm1;
  double dMp, norm_d, norm_dp1, norm_p;
  double obj, objn;

  const int nn = pd.num_free_vertices(); //const int nn = mesh->nn;
  int cg_iter;
  bool valid;

  mHess.initialize( pd, err );  //hMesh(mesh);
  valid = func.update( pd, obj, mGrad, mHess, err ); MSQ_ERRRTN(err);
  if (!valid) {
    MSQ_SETERR(err)("Initial objective function is not valid", MsqError::INVALID_MESH);
    return;
  }
  compute_preconditioner( err ); MSQ_ERRRTN(err);
  pd.recreate_vertices_memento( mMemento, err ); MSQ_ERRRTN(err);

  wVect.resize(nn); w = &wVect[0];//w  = (double *)malloc(3*sizeof(double)*nn);
  zVect.resize(nn); z = &zVect[0];//z  = (double *)malloc(3*sizeof(double)*nn);
  dVect.resize(nn); d = &dVect[0];//d  = (double *)malloc(3*sizeof(double)*nn);
  pVect.resize(nn); p = &pVect[0];//p  = (double *)malloc(3*sizeof(double)*nn);
  rVect.resize(nn); r = &rVect[0];//r  = (double *)malloc(3*sizeof(double)*nn);

  norm_r = norm(&mGrad[0], nn);
  norm_g = sqrt(norm_r);

  while (!term.terminate() && (radius > 1e-20)) {

    memset(d, 0, 3*sizeof(double)*nn);
    memcpy(r, &mGrad[0], nn*sizeof(Vector3D)); //memcpy(r, mesh->g, 3*sizeof(double)*nn);
    norm_g *= cg_tol;

    apply_preconditioner( z, r, err); MSQ_ERRRTN(err); //prec->apply(z, r, prec, mesh);
    negate(p, z, nn);
    rz = inner(r, z, nn);

    dMp    = 0;
    norm_p = rz;
    norm_d = 0;

    cg_iter = 0;
    while ((sqrt(norm_r) > norm_g) && 
#ifdef DO_STEEP_DESC
         (norm_d > tr_num_tol) && 
#endif
         (cg_iter < max_cg_iter)) 
    {
      ++cg_iter;

      memset(w, 0, 3*sizeof(double)*nn);
      matmul(w, mHess, p); //matmul(w, mesh, p);

      kappa = inner(p, w, nn);
      if (kappa <= 0.0) {
        alpha = (sqrt(dMp*dMp+norm_p*(radius-norm_d))-dMp)/norm_p;
        axpy(d, d, alpha, p, nn);
	break;
      }

      alpha = rz / kappa;

      norm_dp1 = norm_d + 2.0*alpha*dMp + alpha*alpha*norm_p;
      if (norm_dp1 >= radius) {
        alpha = (sqrt(dMp*dMp+norm_p*(radius-norm_d))-dMp)/norm_p;
        axpy(d, d, alpha, p, nn);
	break;
      }

      axpy(d, d, alpha, p, nn);
      axpy(r, r, alpha, w, nn);
      norm_r = norm(r, nn);

      apply_preconditioner( z, r, err); MSQ_ERRRTN(err); //prec->apply(z, r, prec, mesh);

      rzm1 = rz;
      rz = inner(r, z, nn);
      beta = rz / rzm1;
      maxpy(p, z, beta, p, nn);

      dMp = beta*(dMp + alpha*norm_p);
      norm_p = rz + beta*beta*norm_p;
      norm_d = norm_dp1;
    }

#ifdef DO_STEEP_DESC    
    if (norm_d <= tr_num_tol) {
      norm_g = norm(&mGrad[0], nn);
      double ll = 1.0;
      if (norm_g < tr_num_tol)
        break;
      if (norm_g > radius*radius)
        ll = radius / sqrt(norm_g);
      for (int i = 0; i < nn; ++i)
        d[i] = ll * mGrad[i];
    }
#endif

    alpha = inner( &mGrad[0], d, nn ); // inner(mesh->g, d, nn);

    memset(p, 0, 3*sizeof(double)*nn);
    matmul(p, mHess, d); //matmul(p, mesh, d);
    beta = 0.5*inner(p, d, nn);
    kappa = alpha + beta;

    /* Put the new point into the locations */
    pd.move_free_vertices_constrained( d, nn, 1.0, err ); MSQ_ERRRTN(err);

    valid = func.evaluate( pd, objn, err ); MSQ_ERRRTN(err); 
    if (!valid) {
      /* Function not defined at trial point */
      radius *= tr_decr_undef;
      pd.set_to_vertices_memento( mMemento, err ); MSQ_ERRRTN(err); 
      continue;
    }
      

    if ((fabs(kappa) <= tr_num_tol) && (fabs(objn - obj) <= tr_num_tol)) {
      kappa = 1;
    }
    else {
      kappa = (objn - obj) / kappa;
    }
    
    if (kappa < eta_1) {
      /* Iterate is unacceptable */
      radius *= tr_decr_def;
      pd.set_to_vertices_memento( mMemento, err ); MSQ_ERRRTN(err); 
      continue;
    }

      /* Iterate is acceptable */
    if (kappa >= eta_2) {
      /* Iterate is a very good step, increase radius */
      radius *= tr_incr;
      if (radius > 1e20) {
	radius = 1e20;
      }
    }

    func.update( pd, obj, mGrad, mHess, err );
    compute_preconditioner( err ); MSQ_ERRRTN(err);
    pd.recreate_vertices_memento( mMemento, err ); MSQ_ERRRTN(err);

    norm_r = norm(&mGrad[0], nn);
    norm_g = sqrt(norm_r);

    // checks stopping criterion 
    term.accumulate_patch( pd, err ); MSQ_ERRRTN(err);
    term.accumulate_inner( pd, objn, &mGrad[0], err ); MSQ_ERRRTN(err);
  }
}    

} // namespace Mesquite
