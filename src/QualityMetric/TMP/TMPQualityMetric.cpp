/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TMPQualityMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#undef PRINT_INFO

#include "Mesquite.hpp"
#include "TMPQualityMetric.hpp"
#include "MsqMatrix.hpp"
#include "ElementQM.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MappingFunction.hpp"
#include "WeightCalculator.hpp"
#include "TargetCalculator.hpp"
#include "TargetMetric2D.hpp"
#include "TargetMetric3D.hpp"
#include "TargetMetricUtil.hpp"

#ifdef PRINT_INFO
#  include <iostream>
#endif

#include <functional>
#include <algorithm>

namespace MESQUITE_NS {


#ifdef PRINT_INFO
template <int R, int C>
void write_vect( char n, const MsqMatrix<R,C>& M )
{
  std::cout << "  " << n << ':';
  for (int c = 0; c < C; ++c) {
    std::cout << '[';
    for (int r = 0; r < R; ++r)
      std::cout << M(r,c) << ' ';
    std::cout << ']';
  }
  std::cout << std::endl;
}

template <int D>
void print_info( size_t elem, Sample sample,
                 const MsqMatrix<3,D>& A,
                 const MsqMatrix<3,D>& W,
                 const MsqMatrix<D,D>& T )
{
  std::cout << "Elem " << elem << " Dim " << sample.dimension << " Num " << sample.number << " :" << std::endl;
  write_vect<3,D>( 'A', A );
  write_vect<3,D>( 'W', W );
  write_vect<D,D>( 'T', T );
}
#endif

int TMPQualityMetric::get_negate_flag( ) const { return 1; }

std::string TMPQualityMetric::get_name() const
{
  std::string result( "TMP(" );
  std::string pfx;
  if (weightCalc) {
    pfx = "c_k ";
  }
  
  if (metric2D && metric3D) {
    std::string name2d = metric2D->get_name();
    std::string name3d = metric3D->get_name();
    if (name2d == name3d)
      result += pfx + name2d;
    else
      result += std::string("2D:") + pfx + name2d + ",3D:" + pfx + name3d;
  }
  else if (metric2D)
    result += pfx + metric2D->get_name();
  else if (metric3D)
    result += pfx + metric3D->get_name();
  else
    result += "NULL";
  
  result += ')';
  return result;
}


void TMPQualityMetric::get_evaluations( PatchData& pd,
                                      std::vector<size_t>& handles,
                                      bool free,
                                      MsqError& err )
{
  get_sample_pt_evaluations( pd, handles, free, err );
}

void TMPQualityMetric::get_patch_evaluations( PatchData& pd,
                                      std::vector<size_t>& handles,
                                      bool free,
                                      MsqError& err )
{
  get_sample_pt_evaluations( pd, handles, free, err );
}

void TMPQualityMetric::get_element_evaluations( PatchData& pd,
                                              size_t elem,
                                              std::vector<size_t>& handles,
                                              MsqError& err )
{
  get_elem_sample_points( pd, elem, handles, err );
}

/**\brief Calculate gradient from derivatives of mapping function terms
 *        and derivatives of target metric. */
template <int DIM> inline
void gradient( size_t num_free_verts,
               const MsqVector<DIM>* dNdxi,
               const MsqMatrix<3,DIM>& dmdA,
               std::vector<Vector3D>& grad )
{
  grad.clear();
  grad.resize( num_free_verts, Vector3D(0,0,0) );
  for (size_t i = 0; i < num_free_verts; ++i)
    grad[i] = Vector3D( (dmdA * dNdxi[i]).data() );
}

/**\brief Calculate Hessian from derivatives of mapping function terms
 *        and derivatives of target metric. */
template <int DIM,typename MAT> inline
void hessian( size_t num_free_verts,
              const MsqVector<DIM>* dNdxi,
              const MsqMatrix<DIM,DIM>* d2mdA2,
              MAT* hess )
{
  MsqMatrix<1,DIM> tmp[DIM][DIM];
  size_t h = 0; // index of current Hessian block

  for (size_t i = 0; i < num_free_verts; ++i) {
  
      // Populate TMP with vector-matrix procucts common
      // to terms of this Hessian row.
    const MsqMatrix<1,DIM>& gi = transpose(dNdxi[i]);
    switch (DIM) {
      case 3:
        tmp[0][2] = gi * d2mdA2[2];
        tmp[1][2] = gi * d2mdA2[4];
        tmp[2][0] = gi * transpose(d2mdA2[2]);
        tmp[2][1] = gi * transpose(d2mdA2[4]);
        tmp[2][2] = gi * d2mdA2[5];
     case 2:
        tmp[0][1] = gi * d2mdA2[1];
        tmp[1][0] = gi * transpose(d2mdA2[1]);
        tmp[1][1] = gi * d2mdA2[DIM];
      case 1:
        tmp[0][0] = gi * d2mdA2[0];
      case 0: 
        break;
      default: assert(false);
    }

      // Calculate Hessian diagonal block
    MAT& H = hess[h++];
    switch (DIM) {
      case 3:
        H(0,2) = H(2,0) = tmp[0][2] * transpose(gi);
        H(1,2) = H(2,1) = tmp[1][2] * transpose(gi);
        H(2,2) =          tmp[2][2] * transpose(gi);
      case 2:
        H(0,1) = H(1,0) = tmp[0][1] * transpose(gi);
        H(1,1) =          tmp[1][1] * transpose(gi);
      case 1:
        H(0,0) =          tmp[0][0] * transpose(gi);
      case 0: 
        break;
      default: assert(false);
    }
    
      // Calculate remainder of Hessian row
    for (size_t j = i+1; j < num_free_verts; ++j) {
      MAT& H = hess[h++];
      const MsqMatrix<DIM,1>& gj = dNdxi[j];
      switch (DIM) {
        case 3:
          H(0,2) = tmp[0][2] * gj;
          H(1,2) = tmp[1][2] * gj;
          H(2,0) = tmp[2][0] * gj;
          H(2,1) = tmp[2][1] * gj;
          H(2,2) = tmp[2][2] * gj;
        case 2:
          H(0,1) = tmp[0][1] * gj;
          H(1,0) = tmp[1][0] * gj;
          H(1,1) = tmp[1][1] * gj;
        case 1:
          H(0,0) = tmp[0][0] * gj;
        case 0: 
          break;
        default: assert(false);
      }
    }
  }
}

bool TMPQualityMetric::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  size_t num_idx;
  return evaluate_with_indices( pd, handle, value, mIndices, num_idx, err );
}


bool TMPQualityMetric::evaluate_with_indices( PatchData& pd,
                                              size_t handle,
                                              double& value,
                                              std::vector<size_t>& indices,
                                              MsqError& err )
{
  indices.resize( MAX_ELEM_NODES );
  size_t num_idx = 0;
  bool result = evaluate_with_indices( pd, handle, value, &indices[0], num_idx, err );
  indices.resize( num_idx );
  return result;
}
/*
static void get_u_perp( const MsqVector<3>& u,
                        MsqVector<3>& u_perp )
{
  double a = sqrt(u[0]*u[0] + u[1]*u[1]);
  if (a < 1e-10) {
    u_perp[0] = 1.0;
    u_perp[1] = u_perp[2] = 0.0;
  }
  else {
    double b = -u[2]/a;
    u_perp[0] = u[0]*b;
    u_perp[1] = u[1]*b;
    u_perp[2] = a;
  }
}
*/
/* Do transform M_hat = S_a M_{3x2}, M_{2x2} Theta^-1 M_hat
 * where the plane into which we are projecting is the cross
 * product of the columns of M, such that S_a is I.  Use the
 * first column of M as u_perp.  
 *
 * Also pass back the cross product of the columns of M as u,
 * and the first column of M as u_perp, both normalized.
 */
static inline void
project_to_matrix_plane( const MsqMatrix<3,2>& M_in,
                         MsqMatrix<2,2>& M_out,
                         MsqVector<3>& u,
                         MsqVector<3>& u_perp )
{
  u = M_in.column(0) * M_in.column(1);
  double u_len = length(u);
  u *= 1.0/u_len;
  u_perp = M_in.column(0);
  double len0 = length(u_perp);
  u_perp *= 1.0/len0;

   // M_out = transpose(theta)*M_in
  M_out(0,0) = len0;
  M_out(0,1) = u_perp % M_in.column(1);
  M_out(1,0) = 0.0;
  M_out(1,1) = u_len / len0;
}

/* Do transform M_hat = S_a M_{3x2}, M_{2x2} Theta^-1 M_hat
 * where the plane into which we are projecting is orthogonal
 * to the passed u vector.
 */
static inline void
project_to_perp_plane(  MsqMatrix<3,2> J,
                        const MsqVector<3>& u,
                        const MsqVector<3>& u_perp,
                        MsqMatrix<2,2>& A,
                        MsqMatrix<3,2>& S_a_transpose_Theta )
{
  MsqVector<3> n_a = J.column(0) * J.column(1);
  n_a *= 1.0/length(n_a);
  double ndot = n_a % u;
  double sigma = (ndot < 0.0) ? -1 : 1;
  double cosphi = sigma * ndot;
  MsqVector<3> cross = n_a * u;
  double sinphi = length(cross);

  MsqMatrix<3,2> Theta;
  Theta.set_column(0,   u_perp);
  Theta.set_column(1, u*u_perp);

    // If columns of J are not in plane orthogonal to u, then
    // rotate J such that they are.
  if (sinphi > 1e-12) {
    MsqVector<3> m = sigma * cross;
    MsqVector<3> n = (1/sinphi) * m;
    MsqVector<3> p = (1-cosphi) * n;
    double s_a[] = 
      { p[0]*n[0] + cosphi, p[0]*n[1] - m[2],   p[0]*n[2] + m[1],
        p[1]*n[0] + m[2],   p[1]*n[1] + cosphi, p[1]*n[2] - m[0],
        p[2]*n[0] - m[1],   p[2]*n[1] + m[0],   p[2]*n[2] + cosphi };
    MsqMatrix<3,3> S_a(s_a);
    J = S_a * J;
    S_a_transpose_Theta = transpose(S_a) * Theta;
  } 
  else {
    S_a_transpose_Theta = Theta;
//    J *= sigma;
  }

    // Project to get 2x2 A from A_hat (which might be equal to J)
  A = transpose(Theta) * J;
}

inline void
TMPQualityMetric::evaluate_surface_common( PatchData& pd,
                                           Sample s,
                                           size_t e,
                                           const NodeSet& bits,
                                           size_t* indices,
                                           size_t& num_indices,
                                           MsqVector<2>* derivs,
                                           MsqMatrix<2,2>& W,
                                           MsqMatrix<2,2>& A,
                                           MsqMatrix<3,2>& S_a_transpose_Theta,
                                           MsqError& err )
{
  EntityTopology type = pd.element_by_index( e ).get_element_type();

  if (!metric2D) {
    MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  const MappingFunction2D* mf = pd.get_mapping_function_2D( type );
  if (!mf) {
    MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }

  MsqMatrix<3,2> J;
  mf->jacobian( pd, e, bits, s, indices, derivs, num_indices, J, err );

    // If we have a 3x2 target matrix 
  if (targetCalc->have_surface_orient()) {
    MsqVector<3> u, u_perp;
    MsqMatrix<3,2> W_hat;
    targetCalc->get_surface_target( pd, e, s, W_hat, err ); MSQ_ERRRTN(err);
      // Use the cross product of the columns of W as the normal of the 
      // plane to work in (i.e. u.).  W should have been constructed such
      // that said cross product is in the direction of (n_s)_init.  And if
      // for some reason it as not, then using something other than said
      // cross product is likely to produce very wrong results.
    project_to_matrix_plane( W_hat, W, u, u_perp );
      // Do the transforms on A to align it with W and project into the plane.
    project_to_perp_plane( J, u, u_perp, A, S_a_transpose_Theta );
  }
    // Otherwise if we have a 2x2 target matrix (i.e. the target does
    // not contain orientation information), project into the plane
    // tangent to J.
  else {
    MsqVector<3> u, u_perp;
    targetCalc->get_2D_target( pd, e, s, W, err ); MSQ_ERRRTN(err);
    project_to_matrix_plane( J, A, u, u_perp );
    S_a_transpose_Theta.set_column(0, u_perp);
    S_a_transpose_Theta.set_column(1, u*u_perp);
      // If the domain is set, adjust the sign of things correctly
      // for the case where the element is inverted with respect
      // to the domain.
    if (pd.domain_set()) {
      Vector3D n;
      pd.get_domain_normal_at_sample( e, s, n, err );
      MSQ_ERRRTN(err);
        // if sigma == -1
      if (Vector3D(u.data()) % n < 0.0) {
          // flip u
        u = -u;
          // S_a_transpose_Theta == Theta, because S_a == I here.
          // u_perp is unaffected by flipping u, so only the second
          // column of S_a_transpose_Theta and the second row of A
          // are flipped because u x u_perp will be flipped.
        S_a_transpose_Theta.set_column(1, -S_a_transpose_Theta.column(1) );
        A.set_row( 1, -A.row(1) );
      }
    }
  }
}                    
                         

bool TMPQualityMetric::evaluate_with_indices( PatchData& pd,
                                              size_t handle,
                                              double& value,
                                              size_t* indices,
                                              size_t& num_indices,
                                              MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const Sample s = ElemSampleQM::sample( handle );
  const size_t e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  unsigned edim = TopologyInfo::dimension( type );
  const NodeSet bits = pd.non_slave_node_set( e );
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    const MappingFunction3D* mf = pd.get_mapping_function_3D( type );
    if (!mf) {
      MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }

    MsqMatrix<3,3> A, W;
    mf->jacobian( pd, e, bits, s, indices, mDerivs3D, num_indices, A, err );
    MSQ_ERRZERO(err);
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    MsqMatrix<2,2> W, A;
    MsqMatrix<3,2> S_a_transpose_Theta;
    evaluate_surface_common( pd, s, e, bits, indices, num_indices,
                             mDerivs2D, W, A, S_a_transpose_Theta, err ); 
                             MSQ_ERRZERO(err);
    rval = metric2D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
#ifdef PRINT_INFO
    print_info<2>( e, s, J, Wp, A * inverse(W) );
#endif
  }
  else {
    assert(false);
    return false;
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }
  
  return rval;
}

bool TMPQualityMetric::evaluate_with_gradient( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           std::vector<size_t>& indices,
                                           std::vector<Vector3D>& grad,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const Sample s = ElemSampleQM::sample( handle );
  const size_t e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  unsigned edim = TopologyInfo::dimension( type );
  size_t num_idx = 0;
  const NodeSet bits = pd.non_slave_node_set( e );
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    const MappingFunction3D* mf = pd.get_mapping_function_3D( type );
    if (!mf) {
      MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }

    MsqMatrix<3,3> A, W, dmdA;
    mf->jacobian( pd, e, bits, s, mIndices, mDerivs3D, num_idx, A, err );
    MSQ_ERRZERO(err);
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_grad( A, W, value, dmdA, err ); MSQ_ERRZERO(err);
    gradient<3>( num_idx, mDerivs3D, dmdA, grad );
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    MsqMatrix<2,2> W, A, dmdA;
    MsqMatrix<3,2> S_a_transpose_Theta;
    evaluate_surface_common( pd, s, e, bits, mIndices, num_idx,
                             mDerivs2D, W, A, S_a_transpose_Theta, err ); 
                             MSQ_ERRZERO(err);
    rval = metric2D->evaluate_with_grad( A, W, value, dmdA, err ); MSQ_ERRZERO(err);
    gradient<2>( num_idx, mDerivs2D, S_a_transpose_Theta*dmdA, grad );
#ifdef PRINT_INFO
    print_info<2>( e, s, J, Wp, A * inverse(W) );
#endif
  }
  else {
    assert(false);
    return false;
  }
  
    // pass back index list
  indices.resize( num_idx );
  std::copy( mIndices, mIndices+num_idx, indices.begin() );
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
    for (size_t i = 0; i < num_idx; ++i)
      grad[i] *= ck;
  }
  
  return rval;
}


bool TMPQualityMetric::evaluate_with_Hessian( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           std::vector<size_t>& indices,
                                           std::vector<Vector3D>& grad,
                                           std::vector<Matrix3D>& Hessian,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const Sample s = ElemSampleQM::sample( handle );
  const size_t e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  unsigned edim = TopologyInfo::dimension( type );
  size_t num_idx = 0;
  const NodeSet bits = pd.non_slave_node_set( e );
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    const MappingFunction3D* mf = pd.get_mapping_function_3D( type );
    if (!mf) {
      MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }

    MsqMatrix<3,3> A, W, dmdA, d2mdA2[6];
    mf->jacobian( pd, e, bits, s, mIndices, mDerivs3D, num_idx, A, err );
    MSQ_ERRZERO(err);
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    gradient<3>( num_idx, mDerivs3D, dmdA, grad );
    Hessian.resize( num_idx*(num_idx+1)/2 );
    if (num_idx)
      hessian<3>( num_idx, mDerivs3D, d2mdA2, &Hessian[0] );
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    MsqMatrix<2,2> W, A, dmdA, d2mdA2[3];
    MsqMatrix<3,2> SaT_Th;
    evaluate_surface_common( pd, s, e, bits, mIndices, num_idx,
                             mDerivs2D, W, A, SaT_Th, err ); 
                             MSQ_ERRZERO(err);
    rval = metric2D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    gradient<2>( num_idx, mDerivs2D, SaT_Th * dmdA, grad );
    const size_t n = num_idx*(num_idx+1)/2;
      // calculate 2D hessian
    hess2d.resize(n);
    if (n)
      hessian<2>( num_idx, mDerivs2D, d2mdA2, &hess2d[0] );
      // calculate surface hessian as transform of 2D hessian
    Hessian.resize(n);
    for (size_t i = 0; i < n; ++i)
      Hessian[i] = Matrix3D( (SaT_Th * hess2d[i] * transpose(SaT_Th)).data() );
#ifdef PRINT_INFO
    print_info<2>( e, s, J, Wp, A * inverse(W) );
#endif
  }
  else {
    assert(0);
    return false;
  }
  
    // pass back index list
  indices.resize( num_idx );
  std::copy( mIndices, mIndices+num_idx, indices.begin() );
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
    for (size_t i = 0; i < num_idx; ++i)
      grad[i] *= ck;
    for (size_t i = 0; i < Hessian.size(); ++i)
      Hessian[i] *= ck;
  }
  
  return rval;
}


bool TMPQualityMetric::evaluate_with_Hessian_diagonal( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           std::vector<size_t>& indices,
                                           std::vector<Vector3D>& grad,
                                           std::vector<SymMatrix3D>& diagonal,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const Sample s = ElemSampleQM::sample( handle );
  const size_t e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  unsigned edim = TopologyInfo::dimension( type );
  size_t num_idx = 0;
  const NodeSet bits = pd.non_slave_node_set( e );
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    const MappingFunction3D* mf = pd.get_mapping_function_3D( type );
    if (!mf) {
      MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }

    MsqMatrix<3,3> A, W, dmdA, d2mdA2[6];
    mf->jacobian( pd, e, bits, s, mIndices, mDerivs3D, num_idx, A, err );
    MSQ_ERRZERO(err);
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    gradient<3>( num_idx, mDerivs3D, dmdA, grad );
    
    diagonal.resize( num_idx );
    for (size_t i = 0; i < num_idx; ++i) {
      SymMatrix3D& H = diagonal[i];
      for (unsigned j = 0; j < 6; ++j)
        H[j] = transpose(mDerivs3D[i]) * d2mdA2[j] * mDerivs3D[i];
    }
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    MsqMatrix<2,2> W, A, dmdA, d2mdA2[3];
    MsqMatrix<3,2> SaT_Th;
    evaluate_surface_common( pd, s, e, bits, mIndices, num_idx,
                             mDerivs2D, W, A, SaT_Th, err ); 
                             MSQ_ERRZERO(err);
    rval = metric2D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    gradient<2>( num_idx, mDerivs2D, SaT_Th * dmdA, grad );

    diagonal.resize( num_idx );
    for (size_t i = 0; i < num_idx; ++i) {
      MsqMatrix<2,2> block2d;
      block2d(0,0) = transpose(mDerivs2D[i]) * d2mdA2[0] * mDerivs2D[i];
      block2d(0,1) = transpose(mDerivs2D[i]) * d2mdA2[1] * mDerivs2D[i];
      block2d(1,0) = block2d(0,1);
      block2d(1,1) = transpose(mDerivs2D[i]) * d2mdA2[2] * mDerivs2D[i];
      MsqMatrix<3,2> p = SaT_Th * block2d;
      
      SymMatrix3D& H = diagonal[i];
      H[0] = p.row(0) * transpose(SaT_Th.row(0));
      H[1] = p.row(0) * transpose(SaT_Th.row(1));
      H[2] = p.row(0) * transpose(SaT_Th.row(2));
      H[3] = p.row(1) * transpose(SaT_Th.row(1));
      H[4] = p.row(1) * transpose(SaT_Th.row(2));
      H[5] = p.row(2) * transpose(SaT_Th.row(2));
    }
#ifdef PRINT_INFO
    print_info<2>( e, s, J, Wp, A * inverse(W) );
#endif
  }
  else {
    assert(0);
    return false;
  }
  
    // pass back index list
  indices.resize( num_idx );
  std::copy( mIndices, mIndices+num_idx, indices.begin() );
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
    for (size_t i = 0; i < num_idx; ++i) {
      grad[i] *= ck;
      diagonal[i] *= ck;
    }
  }
  
  return rval;
}




} // namespace Mesquite
