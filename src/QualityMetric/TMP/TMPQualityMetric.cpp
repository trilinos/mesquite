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

#include "Mesquite.hpp"
#include "TMPQualityMetric.hpp"
#include "SamplePoints.hpp"
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

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <functional.h>
# include <algorithm.h>
#else
# include <functional>
# include <algorithm>
#endif

namespace Mesquite {

int TMPQualityMetric::get_negate_flag( ) const { return 1; }

msq_std::string TMPQualityMetric::get_name() const
  { return msq_std::string("TMPQualityMetric"); }

void TMPQualityMetric::get_evaluations( PatchData& pd,
                                      msq_std::vector<size_t>& handles,
                                      bool free,
                                      MsqError& err )
{
  get_sample_pt_evaluations( pd, samplePts, handles, free, err );
}

void TMPQualityMetric::get_element_evaluations( PatchData& pd,
                                              size_t elem,
                                              msq_std::vector<size_t>& handles,
                                              MsqError& err )
{
  get_elem_sample_points( pd, samplePts, elem, handles, err );
}

void TMPQualityMetric::mapping_function_derivs( PatchData& pd,
                                                size_t handle,
                                                std::vector<size_t>& indices,
                                                std::vector<double>& derivs,
                                                MsqError& err )
{
  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  const EntityTopology type = pd.element_by_index( e ).get_element_type();
  const unsigned edim = TopologyInfo::dimension( type );
  const unsigned ho_bits = pd.higher_order_node_bits( e );
  
  unsigned dim, num;
  samplePts->location_from_sample_number( type, s, dim, num );
 
  const MappingFunction* func = pd.get_mapping_function( type );
  if (!func) {
    MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  indices.clear();
  derivs.clear();
  switch (dim) {
    case 0:
      func->derivatives_at_corner( num, ho_bits, indices, derivs, err );
      break;
    case 1:
      func->derivatives_at_mid_edge( num, ho_bits, indices, derivs, err );
      break;
    case 2:
      if (edim != 2) {
        func->derivatives_at_mid_face( num, ho_bits, indices, derivs, err );
        break;
      }
    case 3:
      func->derivatives_at_mid_elem( ho_bits, indices, derivs, err );
      break;
  }
  MSQ_CHKERR(err);  
}

/** \brief Common code for metric evaluation.  Templatized on dimension so
 *         that it can be used for both 2D and 3D metrics.
 *
 * Given element index and the results of 'mapping_function_derivs'
 * (indices and derivs), do the following:
 * - Calculate Jacobian
 * - Remove entries from indices and derivs corresponding to fixed vertices
 * - As input, 'indices' contains indices into the element connectivity
 *   list.  Change these to indices into the PatchData vertex list.
 */
template <int DIM> inline
void jacobian( PatchData& pd,
               size_t elem_idx,
               std::vector<size_t>& indices,
               std::vector<double>& derivs,
               MsqMatrix<3,DIM>& J )
{
  MsqMeshEntity& elem = pd.element_by_index( elem_idx );
  const size_t* conn = elem.get_vertex_index_array();

  size_t w = 0;
  Vector3D cols[DIM];
  for (size_t i = 0; i < indices.size(); ++i) {
    size_t idx = conn[indices[i]];
    Vector3D coords = pd.vertex_by_index( idx );
    switch (DIM) {
      case 3: cols[2] += derivs[DIM*i+2] * coords;
      case 2: cols[1] += derivs[DIM*i+1] * coords;
      case 1: cols[0] += derivs[DIM*i  ] * coords;
      case 0: break;
      default: assert(false);
    }
    
    if (idx < pd.num_free_vertices()) {
      indices[w] = idx;
      switch (DIM) {
        case 3: derivs[DIM*w+2] = derivs[DIM*i+2];
        case 2: derivs[DIM*w+1] = derivs[DIM*i+1];
        case 1: derivs[DIM*w  ] = derivs[DIM*i  ];
        case 0: break;
      }
      ++w;
    }
  }
  indices.resize(w);
  derivs.resize(DIM*w);
  J = MsqMatrix<3,DIM>( reinterpret_cast< MsqMatrix<3,1>* >(cols) );
}

/**\brief Calculate gradient from derivatives of mapping function terms
 *        and derivatives of target metric. */
template <int DIM> inline
void gradient( size_t num_free_verts,
               const MsqMatrix<DIM,1>* dNdxi,
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
              const MsqMatrix<DIM,1>* dNdxi,
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
  mIndices.clear();
  return evaluate_with_indices( pd, handle, value, mIndices, err );
}


bool TMPQualityMetric::evaluate_with_indices( PatchData& pd,
                                              size_t handle,
                                              double& value,
                                              msq_std::vector<size_t>& indices,
                                              MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  unsigned edim = TopologyInfo::dimension( elem.get_element_type() );
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }

    MsqMatrix<3,3> A, W;
    jacobian<3>( pd, e, indices, mDerivs, A );
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    
    MsqMatrix<3,2> J, Wp, RZ;
    jacobian<2>( pd, e, indices, mDerivs, J );
    targetCalc->get_2D_target( pd, e, samplePts, s, Wp, err ); MSQ_ERRZERO(err);
    
    MsqMatrix<2,2> W, A;
    surface_to_2d( J, Wp, W, RZ );
    A = transpose(RZ) * J;
    rval = metric2D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
  }
  else {
    assert(false);
    return false;
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, samplePts, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }
  
  return rval;
}
                 

bool TMPQualityMetric::evaluate_with_gradient( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           msq_std::vector<size_t>& indices,
                                           msq_std::vector<Vector3D>& grad,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  unsigned edim = TopologyInfo::dimension( elem.get_element_type() );
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
  bool rval;
  std::vector<double>::const_iterator d = mDerivs.begin();
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    MsqMatrix<3,3> A, W, dmdA;
    jacobian<3>( pd, e, indices, mDerivs, A );
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_grad( A, W, value, dmdA, err ); MSQ_ERRZERO(err);
    gradient<3>( indices.size(), 
                 reinterpret_cast< MsqMatrix<3,1>* >(&mDerivs[0]),
                 dmdA, grad );
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    
    MsqMatrix<3,2> J, Wp, RZ;
    jacobian<2>( pd, e, indices, mDerivs, J );
    targetCalc->get_2D_target( pd, e, samplePts, s, Wp, err ); MSQ_ERRZERO(err);
    
    MsqMatrix<2,2> W, dmdA, A;
    surface_to_2d( J, Wp, W, RZ );
    A = transpose(RZ) * J;
    rval = metric2D->evaluate_with_grad( A, W, value, dmdA, err ); MSQ_ERRZERO(err);
    gradient<2>( indices.size(), 
                 reinterpret_cast< MsqMatrix<2,1>* >(&mDerivs[0]),
                 RZ * dmdA, grad );
  }
  else {
    assert(false);
    return false;
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, samplePts, s, err ); MSQ_ERRZERO(err);
    value *= ck;
    for (size_t i = 0; i < grad.size(); ++i)
      grad[i] *= ck;
  }
  
  return rval;
}


bool TMPQualityMetric::evaluate_with_Hessian( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           msq_std::vector<size_t>& indices,
                                           msq_std::vector<Vector3D>& grad,
                                           msq_std::vector<Matrix3D>& Hessian,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  unsigned edim = TopologyInfo::dimension( elem.get_element_type() );
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    MsqMatrix<3,3> A, W, dmdA, d2mdA2[6];
    jacobian<3>( pd, e, indices, mDerivs, A );
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    const MsqMatrix<3,1>* dNdxi = reinterpret_cast< MsqMatrix<3,1>* >(&mDerivs[0]);
    gradient<3>( indices.size(), dNdxi, dmdA, grad );
    Hessian.resize( indices.size()*(indices.size()+1)/2 );
    hessian<3>( indices.size(), dNdxi, d2mdA2, &Hessian[0] );
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    
    MsqMatrix<3,2> J, Wp, RZ;
    jacobian<2>( pd, e, indices, mDerivs, J );
    targetCalc->get_2D_target( pd, e, samplePts, s, Wp, err ); MSQ_ERRZERO(err);
    
    MsqMatrix<2,2> W, dmdA, A, d2mdA2[3];
    surface_to_2d( J, Wp, W, RZ );
    A = transpose(RZ) * J;
    rval = metric2D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    const MsqMatrix<2,1>* dNdxi = reinterpret_cast< MsqMatrix<2,1>* >(&mDerivs[0]);
    gradient<2>( indices.size(), dNdxi, RZ * dmdA, grad );
    const size_t n = indices.size()*(indices.size()+1)/2;
      // calculate 2D hessian
    hess2d.resize(n);
    hessian<2>( indices.size(), dNdxi, d2mdA2, &hess2d[0] );
      // calculate surface hessian as transform of 2D hessian
    Hessian.resize(n);
    for (size_t i = 0; i < n; ++i)
      Hessian[i] = Matrix3D( (RZ * transpose(hess2d[i]) * transpose(RZ)).data() );
  }
  else {
    assert(0);
    return false;
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, samplePts, s, err ); MSQ_ERRZERO(err);
    value *= ck;
    for (size_t i = 0; i < grad.size(); ++i)
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
                                           msq_std::vector<size_t>& indices,
                                           msq_std::vector<Vector3D>& grad,
                                           msq_std::vector<SymMatrix3D>& diagonal,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  unsigned edim = TopologyInfo::dimension( elem.get_element_type() );
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    MsqMatrix<3,3> A, W, dmdA, d2mdA2[6];
    jacobian<3>( pd, e, indices, mDerivs, A );
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    const MsqMatrix<3,1>* dNdxi = reinterpret_cast< MsqMatrix<3,1>* >(&mDerivs[0]);
    gradient<3>( indices.size(), dNdxi, dmdA, grad );
    
    diagonal.resize( indices.size() );
    for (size_t i = 0; i < indices.size(); ++i) {
      SymMatrix3D& H = diagonal[i];
      for (unsigned j = 0; j < 6; ++j)
        H[j] = transpose( dNdxi[i] ) * d2mdA2[j] * dNdxi[i]; 
    }
  }
  else if (edim == 2) {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    
    MsqMatrix<3,2> J, Wp, RZ;
    jacobian<2>( pd, e, indices, mDerivs, J );
    targetCalc->get_2D_target( pd, e, samplePts, s, Wp, err ); MSQ_ERRZERO(err);
    
    MsqMatrix<2,2> W, dmdA, A, d2mdA2[3];
    surface_to_2d( J, Wp, W, RZ );
    A = transpose(RZ) * J;
    rval = metric2D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    const MsqMatrix<2,1>* dNdxi = reinterpret_cast< MsqMatrix<2,1>* >(&mDerivs[0]);
    gradient<2>( indices.size(), dNdxi, RZ * dmdA, grad );

    diagonal.resize( indices.size() );
    for (size_t i = 0; i < indices.size(); ++i) {
      MsqMatrix<2,2> block2d;
      block2d(0,0) = transpose(dNdxi[i]) * d2mdA2[0] * dNdxi[i];
      block2d(0,1) = transpose(dNdxi[i]) * d2mdA2[1] * dNdxi[i];
      block2d(1,0) = block2d(0,1);
      block2d(1,1) = transpose(dNdxi[i]) * d2mdA2[2] * dNdxi[i];
      MsqMatrix<3,2> p = RZ * block2d;
      
      SymMatrix3D& H = diagonal[i];
      H[0] = p.row(0) * transpose(RZ.row(0));
      H[1] = p.row(0) * transpose(RZ.row(1));
      H[2] = p.row(0) * transpose(RZ.row(2));
      H[3] = p.row(1) * transpose(RZ.row(1));
      H[4] = p.row(1) * transpose(RZ.row(2));
      H[5] = p.row(2) * transpose(RZ.row(2));
    }
  }
  else {
    assert(0);
    return false;
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, samplePts, s, err ); MSQ_ERRZERO(err);
    value *= ck;
    for (size_t i = 0; i < indices.size(); ++i) {
      grad[i] *= ck;
      diagonal[i] *= ck;
    }
  }
  
  return rval;
}




} // namespace Mesquite
