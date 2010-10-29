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


/** \file TRelQualityMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#undef PRINT_INFO

#include "Mesquite.hpp"
#include "TRelQualityMetric.hpp"
#include "MsqMatrix.hpp"
#include "ElementQM.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MappingFunction.hpp"
#include "WeightCalculator.hpp"
#include "TargetCalculator.hpp"
#include "TRel2DMetric.hpp"
#include "TRel3DMetric.hpp"
#include "TargetMetricUtil.hpp"
#include "TMPDerivs.hpp"

#ifdef PRINT_INFO
#  include <iostream>
#endif

#include <functional>
#include <algorithm>

namespace MESQUITE_NS {

std::string TRelQualityMetric::get_name() const
{
  return make_name( "TRel", metric2D ? metric2D->get_name() : std::string(),
                            metric3D ? metric3D->get_name() : std::string() );
}

bool TRelQualityMetric::evaluate_with_indices( PatchData& pd,
                                               size_t handle,
                                               double& value,
                                               size_t* indices,
                                               size_t& num_indices,
                                               MsqError& err )
{
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
    const MsqMatrix<3,3> Winv = inverse(W);
    const MsqMatrix<3,3> T = A*Winv;
    rval = metric3D->evaluate( T, value, err ); MSQ_ERRZERO(err);
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    MsqMatrix<2,2> W, A;
    MsqMatrix<3,2> S_a_transpose_Theta;
    rval = evaluate_surface_common( pd, s, e, bits, indices, num_indices,
                                 mDerivs2D, W, A, S_a_transpose_Theta, err ); 
    if (MSQ_CHKERR(err) || !rval)
      return false;
    const MsqMatrix<2,2> Winv = inverse(W);
    const MsqMatrix<2,2> T = A * Winv;
    rval = metric2D->evaluate( T, value, err ); MSQ_ERRZERO(err);
#ifdef PRINT_INFO
    print_info<2>( e, s, J, Wp, A * inverse(W) );
#endif
  }
  else {
    assert(false);
    return false;
  }
  
  return rval;
}

bool TRelQualityMetric::evaluate_with_gradient( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           std::vector<size_t>& indices,
                                           std::vector<Vector3D>& grad,
                                           MsqError& err )
{
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

    MsqMatrix<3,3> A, W, dmdT;
    mf->jacobian( pd, e, bits, s, mIndices, mDerivs3D, num_idx, A, err );
    MSQ_ERRZERO(err);
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    const MsqMatrix<3,3> Winv = inverse(W);
    const MsqMatrix<3,3> T = A*Winv;
    rval = metric3D->evaluate_with_grad( T, value, dmdT, err ); MSQ_ERRZERO(err);
    gradient<3>( num_idx, mDerivs3D, dmdT * transpose(Winv), grad );
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    MsqMatrix<2,2> W, A, dmdT;
    MsqMatrix<3,2> S_a_transpose_Theta;
    rval = evaluate_surface_common( pd, s, e, bits, mIndices, num_idx,
                             mDerivs2D, W, A, S_a_transpose_Theta, err ); 
    if (MSQ_CHKERR(err) || !rval)
      return false;
    const MsqMatrix<2,2> Winv = inverse(W);
    const MsqMatrix<2,2> T = A*Winv;
    rval = metric2D->evaluate_with_grad( T, value, dmdT, err ); MSQ_ERRZERO(err);
    gradient<2>( num_idx, mDerivs2D, S_a_transpose_Theta*dmdT*transpose(Winv), grad );
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
  weight( pd, s, e, num_idx, value, &grad[0], 0, 0, err ); MSQ_ERRZERO(err);
  return rval;
}


bool TRelQualityMetric::evaluate_with_Hessian( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           std::vector<size_t>& indices,
                                           std::vector<Vector3D>& grad,
                                           std::vector<Matrix3D>& Hessian,
                                           MsqError& err )
{
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

    MsqMatrix<3,3> A, W, dmdT, d2mdT2[6];
    mf->jacobian( pd, e, bits, s, mIndices, mDerivs3D, num_idx, A, err );
    MSQ_ERRZERO(err);
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    const MsqMatrix<3,3> Winv = inverse(W);
    const MsqMatrix<3,3> T = A*Winv;
    rval = metric3D->evaluate_with_hess( T, value, dmdT, d2mdT2, err ); MSQ_ERRZERO(err);
    gradient<3>( num_idx, mDerivs3D, dmdT*transpose(Winv), grad );
    second_deriv_wrt_product_factor( d2mdT2, Winv );
    Hessian.resize( num_idx*(num_idx+1)/2 );
    if (num_idx)
      hessian<3>( num_idx, mDerivs3D, d2mdT2, &Hessian[0] );
    
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }

    // return finite difference approximation for now

    return QualityMetric::evaluate_with_Hessian( pd, handle,
                                           value, indices, grad, Hessian,
                                           err );
    /*
    MsqMatrix<2,2> W, A, dmdT, d2mdT2[3];
    MsqMatrix<3,2> SaT_Th;
    rval = evaluate_surface_common( pd, s, e, bits, mIndices, num_idx,
                             mDerivs2D, W, A, SaT_Th, err ); 
    if (MSQ_CHKERR(err) || !rval)
      return false;
    const MsqMatrix<2,2> Winv = inverse(W);
    const MsqMatrix<2,2> T = A*Winv;
    rval = metric2D->evaluate_with_hess( T, value, dmdT, d2mdT2, err ); MSQ_ERRZERO(err);
    gradient<2>( num_idx, mDerivs2D, SaT_Th * dmdT * transpose(Winv), grad );
    second_deriv_wrt_product_factor( d2mdT2, Winv );
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
    */
  }
  else {
    assert(0);
    return false;
  }
  
    // pass back index list
  indices.resize( num_idx );
  std::copy( mIndices, mIndices+num_idx, indices.begin() );
  
    // apply target weight to value
  weight( pd, s, e, num_idx, value, &grad[0], 0, &Hessian[0], err ); MSQ_ERRZERO(err);
  return rval;
}


bool TRelQualityMetric::evaluate_with_Hessian_diagonal( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           std::vector<size_t>& indices,
                                           std::vector<Vector3D>& grad,
                                           std::vector<SymMatrix3D>& diagonal,
                                           MsqError& err )
{
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

    MsqMatrix<3,3> A, W, dmdT, d2mdT2[6];
    mf->jacobian( pd, e, bits, s, mIndices, mDerivs3D, num_idx, A, err );
    MSQ_ERRZERO(err);
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    const MsqMatrix<3,3> Winv = inverse(W);
    const MsqMatrix<3,3> T = A*Winv;
    rval = metric3D->evaluate_with_hess( T, value, dmdT, d2mdT2, err ); MSQ_ERRZERO(err);
    gradient<3>( num_idx, mDerivs3D, dmdT * transpose(Winv), grad );
    second_deriv_wrt_product_factor( d2mdT2, Winv );
    
    diagonal.resize( num_idx );
    hessian_diagonal<3>(num_idx, mDerivs3D, d2mdT2, &diagonal[0] );
#ifdef PRINT_INFO
    print_info<3>( e, s, A, W, A * inverse(W) );
#endif
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }

    // use finite diference approximation for now
    return QualityMetric::evaluate_with_Hessian_diagonal( pd, handle,
                                           value, indices, grad, diagonal,
                                           err );
/*
    MsqMatrix<2,2> W, A, dmdT, d2mdT2[3];
    MsqMatrix<3,2> SaT_Th;
    rval = evaluate_surface_common( pd, s, e, bits, mIndices, num_idx,
                             mDerivs2D, W, A, SaT_Th, err ); 
    if (MSQ_CHKERR(err) || !rval)
      return false;
    const MsqMatrix<2,2> Winv = inverse(W);
    const MsqMatrix<2,2> T = A*Winv;
    rval = metric2D->evaluate_with_hess( T, value, dmdT, d2mdT2, err ); MSQ_ERRZERO(err);
    gradient<2>( num_idx, mDerivs2D, SaT_Th * dmdT * transpose(Winv), grad );
    second_deriv_wrt_product_factor( d2mdT2, Winv );

    diagonal.resize( num_idx );
    for (size_t i = 0; i < num_idx; ++i) {
      MsqMatrix<2,2> block2d;
      block2d(0,0) = transpose(mDerivs2D[i]) * d2mdT2[0] * mDerivs2D[i];
      block2d(0,1) = transpose(mDerivs2D[i]) * d2mdT2[1] * mDerivs2D[i];
      block2d(1,0) = block2d(0,1);
      block2d(1,1) = transpose(mDerivs2D[i]) * d2mdT2[2] * mDerivs2D[i];
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
*/
  }
  else {
    assert(0);
    return false;
  }
  
    // pass back index list
  indices.resize( num_idx );
  std::copy( mIndices, mIndices+num_idx, indices.begin() );
  
    // apply target weight to value
  weight( pd, s, e, num_idx, value, &grad[0], &diagonal[0], 0, err ); MSQ_ERRZERO(err);
  return rval;
}



} // namespace Mesquite
