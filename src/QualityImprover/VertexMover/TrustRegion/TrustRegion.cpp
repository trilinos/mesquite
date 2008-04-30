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

  //Precond  *prec;

  //double *v;
  Vector3D *w;
  Vector3D *z;
  Vector3D *d;
  Vector3D *p;
  Vector3D *r;
  //Vector3D *t;

  double norm_r, norm_g;
  double alpha, beta, kappa;
  double rz, rzm1;
  double dMp, norm_d, norm_dp1, norm_p;
  double obj, objn;

  const int nn = pd.num_free_vertices(); //const int nn = mesh->nn;
  int iter = 0, cg_iter;
  int nf = 1, ng = 1, nh = 0;
  int needH = 1;
  int ferr;

  //double m, m1, m2;
  //double s, s1, s2;
  //double t1, t2;

  //struct rusage r0, r1;

#ifdef DISPLAY_MAX
  double fmax;
  int    fidx;
#endif

  if (nn <= 0) {
    /* No nodes!  Just return.                                               */
    return; //return 0;
  }

  //getrusage(RUSAGE_SELF, &r0);

#ifdef USE_WEIGHT
  mesh->w = (double *)malloc(16*sizeof(double)*mesh->ne);

  {
    double *weight = mesh->w;

    double W0[3];
    double W1[3];
    double W2[3];
    double W3[3];
    double Wmat[3][3];
    double Q[3][3];
    double R[3][3];
    double Rinv[3][3];

    int i;

    srand48(1001);

    for (i = 0; i < mesh->ne; ++i) {

#ifdef RANDOM_WEIGHT
      W0[0] = 0.0;
      W0[1] = 0.0;
      W0[2] = 0.0;

      W1[0] = 5.0*drand48() + 0.5;
      W1[1] = 0.0;
      W1[2] = 0.0;

      W2[0] = 10.0*drand48() - 5.0;
      W2[1] = 5.0*drand48() + 0.5;
      W2[2] = 0.0;

      W3[0] = 10.0*drand48() - 5.0;
      W3[1] = 10.0*drand48() - 5.0;
      W3[2] = 5.0*drand48() + 0.5;
#else
      W0[0] = 0.0;
      W0[1] = 0.0;
      W0[2] = 0.0;

      W1[0] = 1.0;
      W1[1] = 0.0;
      W1[2] = 0.0;

      W2[0] = 0.5;
      W2[1] = 0.5 * sqrt(3.0);
      W2[2] = 0.0;

      W3[0] = 0.5;
      W3[1] = 0.5 / sqrt(3.0);
      W3[2] = sqrt(6.0) / 3.0;
#endif

      getmat(Wmat, W0, W1, W2, W3);
      getQR(Q, R, Wmat);
      getinv(Rinv, R);

      weight[0] = (fabs(Rinv[0][0]) > 1E-10) ? Rinv[0][0] : 0.0;
      weight[1] = (fabs(Rinv[0][1]) > 1E-10) ? Rinv[0][1] : 0.0;
      weight[2] = (fabs(Rinv[0][2]) > 1E-10) ? Rinv[0][2] : 0.0;
      weight[3] = (fabs(Rinv[1][1]) > 1E-10) ? Rinv[1][1] : 0.0;
      weight[4] = (fabs(Rinv[1][2]) > 1E-10) ? Rinv[1][2] : 0.0;
      weight[5] = (fabs(Rinv[2][2]) > 1E-10) ? Rinv[2][2] : 0.0;

      weight[6] = (fabs(Q[0][0]) > 1E-10) ? Q[0][0] : 0.0;
      weight[7] = (fabs(Q[0][1]) > 1E-10) ? Q[0][1] : 0.0;
      weight[8] = (fabs(Q[0][2]) > 1E-10) ? Q[0][2] : 0.0;
      weight[9] = (fabs(Q[1][0]) > 1E-10) ? Q[1][0] : 0.0;
      weight[10] = (fabs(Q[1][1]) > 1E-10) ? Q[1][1] : 0.0;
      weight[11] = (fabs(Q[1][2]) > 1E-10) ? Q[1][2] : 0.0;
      weight[12] = (fabs(Q[2][0]) > 1E-10) ? Q[2][0] : 0.0;
      weight[13] = (fabs(Q[2][1]) > 1E-10) ? Q[2][1] : 0.0;
      weight[14] = (fabs(Q[2][2]) > 1E-10) ? Q[2][2] : 0.0;

      weight += 16;
    }
  }
#endif

  mHess.initialize( pd, err );  //hMesh(mesh);
  if (!func.update( pd, obj, mGrad, err ) && MSQ_CHKERR(err)) { //if (gFcn(&obj, mesh)) {
    return ; //  fprintf(stderr, "Invalid starting point.\n");
             //  exit(-1);
  }

  //prec = preCreate(precond, nn, mesh->nz);

  //v  = (double *)malloc(3*sizeof(double)*nn);
  wVect.resize(nn); w = &wVect[0];//w  = (double *)malloc(3*sizeof(double)*nn);
  zVect.resize(nn); z = &zVect[0];//z  = (double *)malloc(3*sizeof(double)*nn);
  dVect.resize(nn); d = &dVect[0];//d  = (double *)malloc(3*sizeof(double)*nn);
  pVect.resize(nn); p = &pVect[0];//p  = (double *)malloc(3*sizeof(double)*nn);
  rVect.resize(nn); r = &rVect[0];//r  = (double *)malloc(3*sizeof(double)*nn);

  //gatherMesh(v, mesh);
  norm_r = norm(&mGrad[0], nn);
  norm_g = sqrt(norm_r);

  //getrusage(RUSAGE_SELF, &r1);

  //m1 = (double) r1.ru_utime.tv_usec;
  //m2 = (double) r0.ru_utime.tv_usec;
  //m = m1 - m2;
    
  //s1 = (double) r1.ru_utime.tv_sec;
  //s2 = (double) r0.ru_utime.tv_sec;
  //s = s1 - s2;

  //t1 = s + m / MICROSEC;

  //m1 = (double) r1.ru_stime.tv_usec;
  //m2 = (double) r0.ru_stime.tv_usec;
  //m = m1 - m2;

  //s1 = (double) r1.ru_stime.tv_sec;
  //s2 = (double) r0.ru_stime.tv_sec;
  //s = s1 - s2;

  //t2 = s + m / MICROSEC;

#ifdef DISPLAY_MAX
  //oMax(&fmax, &fidx, mesh);
  //printf("%4d I      %10.9e %5.4e            %5.4e %5d "
  //       "usr: %5.4e sys: %5.4e tot: %5.4e\n",
  //       iter, obj, norm_g, fmax, fidx, t1, t2, t1+t2);
#else
  //printf("%4d I      %10.9e %5.4e            "
  //       "usr: %5.4e sys: %5.4e tot: %5.4e\n",
  //       iter, obj, norm_g, t1, t2, t1+t2);
#endif

  while (!term.terminate() && (radius > 1e-20)) {
    pd.recreate_vertices_memento( mMemento, err ); MSQ_ERRRTN(err);

    if (needH) {
      //getrusage(RUSAGE_SELF, &r0);
      //++iter;

      //++nh;
      func.update( pd, obj, mGrad, mHess, err ); MSQ_ERRRTN(err); //hOnly(mesh);
      compute_preconditioner( err ); MSQ_ERRRTN(err); //prec->calc(prec, mesh);
      needH = 0;
    }

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
    //axpy(w, v, 1.0, d, nn);
    //scatterMesh(mesh, w);
    pd.move_free_vertices_constrained( d, nn, 1.0, err ); MSQ_ERRRTN(err);

    ++nf;
    ferr = func.evaluate( pd, objn, err ); MSQ_ERRRTN(err); //ferr = oFcn(&objn, mesh);

#ifdef OUTPUT_STATS
    printf("alpha = %20.19e\n", alpha);
    printf("beta  = %20.19e\n", beta);
    printf("kappa = %20.19e\n", kappa);

    if (!ferr) {
      printf("new   = %20.19e\n", objn);
      printf("old   = %20.19e\n", obj);
      printf("ared  = %20.19e\n", objn - obj);
      printf("pred  = %20.19e\n", kappa);
      printf("ratio = %20.19e\n", (objn - obj) / kappa);
    }
#endif

    if (ferr) {
      if ((fabs(kappa) <= tr_num_tol) && (fabs(objn - obj) <= tr_num_tol)) {
	kappa = 1;
      }
      else {
	kappa = (objn - obj) / kappa;
      }

      if (kappa >= eta_1) {
	/* Iterate is acceptable */

        if (kappa >= eta_2) {
	  /* Iterate is a very good step, increase radius */
	  radius *= tr_incr;
	  if (radius > 1e20) {
	    radius = 1e20;
	  }
	}

	/* Update the iterate (v = current point, w = new point) so swap */
	//t = v;
	//v = w;
	//w = t;

	++ng;
	obj = objn;
	func.evaluate( pd, objn, mGrad, err ); //gOnly(mesh);
	needH = 1;
      }
      else {
	/* Iterate is unacceptable */
	radius *= tr_decr_def;
        pd.set_to_vertices_memento( mMemento, err ); MSQ_ERRRTN(err); //scatterMesh(mesh, v);
      }
    }
    else {
      /* Function not defined at trial point */
      radius *= tr_decr_undef;
      pd.set_to_vertices_memento( mMemento, err ); MSQ_ERRRTN(err); //scatterMesh(mesh, v);
    }

    norm_r = norm(&mGrad[0], nn);
    norm_g = sqrt(norm_r);

    //if (needH) {
      //getrusage(RUSAGE_SELF, &r1);

      //m1 = (double) r1.ru_utime.tv_usec;
      //m2 = (double) r0.ru_utime.tv_usec;
      //m = m1 - m2;
    
      //s1 = (double) r1.ru_utime.tv_sec;
      //s2 = (double) r0.ru_utime.tv_sec;
      //s = s1 - s2;
      
      //t1 = s + m / MICROSEC;
      
      //m1 = (double) r1.ru_stime.tv_usec;
      //m2 = (double) r0.ru_stime.tv_usec;
      //m = m1 - m2;
    
      //s1 = (double) r1.ru_stime.tv_sec;
      //s2 = (double) r0.ru_stime.tv_sec;
      //s = s1 - s2;
      
      //t2 = s + m / MICROSEC;
      
#ifdef DISPLAY_MAX
      //oMax(&fmax, &fidx, mesh);
      //printf("%4d N %4d %10.9e %5.4e %5.4e %5.4e %5d "
      //       "usr: %5.4e sys: %5.4e tot: %5.4e\n",
      //       iter, cg_iter, obj, norm_g, radius, fmax, fidx, t1, t2, t1+t2);
#else
      //printf("%4d N %4d %10.9e %5.4e %5.4e usr: %5.4e sys: %5.4e tot: %5.4e\n",
      //       iter, cg_iter, obj, norm_g, radius, t1, t2, t1+t2);
#endif
    //}

    // checks stopping criterion 
    term.accumulate_patch( pd, err ); MSQ_ERRRTN(err);
    term.accumulate_inner( pd, objn, &mGrad[0], err ); MSQ_ERRRTN(err);
  }

  //printf("Function evals: %4d\nGradient evals: %4d\nHessian  evals: %4d\n", 
  //	 nf, ng, nh);

  //prec->destroy(prec);

  //free(v);
  //free(w);
  //free(z);
  //free(d);
  //free(p);
  //free(r);
  //return 0;
}    

} // namespace Mesquite
