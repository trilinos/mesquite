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
  const size_t* conn = elem.get_vertex_index_array();
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
    // Convert from indices into element connectivity list to
    // indices into vertex array in patch data.
  for (msq_std::vector<size_t>::iterator i = indices.begin(); i != indices.end(); ++i)
    *i = conn[*i];
  
  bool rval;
  std::vector<double>::const_iterator d = mDerivs.begin();
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    Vector3D c[3] = { Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0) };
    for (size_t i = 0; i < indices.size(); ++i) {
      Vector3D coords = pd.vertex_by_index( indices[i] );
      c[0] += *d * coords; ++d;
      c[1] += *d * coords; ++d;
      c[2] += *d * coords; ++d;
    }
    MsqMatrix<3,3> A( (MsqMatrix<3,1>*)c );

    MsqMatrix<3,3> W;
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
  }
  else {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    Vector3D c[2] = { Vector3D(0,0,0), Vector3D(0,0,0) };
    for (size_t i = 0; i < indices.size(); ++i) {
      Vector3D coords = pd.vertex_by_index( indices[i] );
      c[0] += *d * coords; ++d;
      c[1] += *d * coords; ++d;
    }
    MsqMatrix<3,2> J( (MsqMatrix<3,1>*)c );
    
    
    MsqMatrix<3,2> Wp;
    targetCalc->get_2D_target( pd, e, samplePts, s, Wp, err ); MSQ_ERRZERO(err);
    
    MsqMatrix<2,2> W;
    MsqMatrix<3,2> RZ;
    surface_to_2d( J, Wp, W, RZ );
    MsqMatrix<2,2> A = transpose(RZ) * J;
    rval = metric2D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, samplePts, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }
  
    // remove indices for non-free vertices
  indices.erase( msq_std::remove_if( indices.begin(), indices.end(), 
    msq_std::bind2nd(msq_std::greater_equal<size_t>(),pd.num_free_vertices())),
    indices.end() );
  
  return rval;
}
                 

bool TMPQualityMetric::evaluate_with_gradient( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           msq_std::vector<size_t>& indices,
                                           msq_std::vector<Vector3D>& gradient,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  unsigned edim = TopologyInfo::dimension( elem.get_element_type() );
  const size_t* conn = elem.get_vertex_index_array();
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
  bool rval;
  std::vector<double>::const_iterator d = mDerivs.begin();
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    Vector3D c[3] = { Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0) };
    for (size_t i = 0; i < indices.size(); ++i) {
        // Convert from indices into element connectivity list to
        // indices into vertex array in patch data.
      indices[i] = conn[indices[i]];
        // calculate Jacobian
      Vector3D coords = pd.vertex_by_index( indices[i] );
      c[0] += *d * coords; ++d;
      c[1] += *d * coords; ++d;
      c[2] += *d * coords; ++d;
    }
    MsqMatrix<3,3> A( (MsqMatrix<3,1>*)c );

    MsqMatrix<3,3> W, dmdA;
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_grad( A, W, value, dmdA, err ); MSQ_ERRZERO(err);
    gradient.clear();
    d = mDerivs.begin();
    for (size_t i =0; i < indices.size(); ++i, d += 3) {
        // if vertex is not fixed, include in gradient
      if (indices[i] < pd.num_free_vertices())
        gradient.push_back(( (dmdA * MsqMatrix<3,1>(&*d)).data() ));
    }
  }
  else {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
  
    Vector3D c[2] = { Vector3D(0,0,0), Vector3D(0,0,0) };
    for (size_t i = 0; i < indices.size(); ++i) {
        // Convert from indices into element connectivity list to
        // indices into vertex array in patch data.
      indices[i] = conn[indices[i]];
        // calculate Jacobian
      Vector3D coords = pd.vertex_by_index( indices[i] );
      c[0] += *d * coords; ++d;
      c[1] += *d * coords; ++d;
    }
    MsqMatrix<3,2> J( (MsqMatrix<3,1>*)c );
    
    
    MsqMatrix<3,2> Wp;
    targetCalc->get_2D_target( pd, e, samplePts, s, Wp, err ); MSQ_ERRZERO(err);
    
    MsqMatrix<2,2> W;
    MsqMatrix<3,2> RZ;
    surface_to_2d( J, Wp, W, RZ );
    MsqMatrix<2,2> dmdA, A = transpose(RZ) * J;
    rval = metric2D->evaluate_with_grad( A, W, value, dmdA, err ); MSQ_ERRZERO(err);
  
    const MsqMatrix<3,2> M = RZ * dmdA;
    gradient.clear();
    d = mDerivs.begin();
    for (size_t i =0; i < indices.size(); ++i, d += 2) {
        // if vertex is not fixed, include in gradient
      if (indices[i] < pd.num_free_vertices())
        gradient.push_back(( (M * MsqMatrix<2,1>(&*d)).data() ));
    }
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, samplePts, s, err ); MSQ_ERRZERO(err);
    value *= ck;
    for (size_t i = 0; i < gradient.size(); ++i)
      gradient[i] *= ck;
  }
  
    // remove indices for non-free vertices
  indices.erase( msq_std::remove_if( indices.begin(), indices.end(), 
    msq_std::bind2nd(msq_std::greater_equal<size_t>(),pd.num_free_vertices())),
    indices.end() );
  
  return rval;
}


bool TMPQualityMetric::evaluate_with_Hessian( 
                                           PatchData& pd,
                                           size_t handle,
                                           double& value,
                                           msq_std::vector<size_t>& indices,
                                           msq_std::vector<Vector3D>& gradient,
                                           msq_std::vector<Matrix3D>& Hessian,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  unsigned edim = TopologyInfo::dimension( elem.get_element_type() );
  const size_t* conn = elem.get_vertex_index_array();
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
  if (edim != 3) // use finite difference approximation for surface elements
    return QualityMetric::evaluate_with_Hessian( pd, handle, value, indices, gradient, Hessian, err );
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    Vector3D c[3] = { Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0) };
    size_t w = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        // Convert from indices into element connectivity list to
        // indices into vertex array in patch data.
      size_t vtx = conn[indices[i]];
        // calculate Jacobian
      Vector3D coords = pd.vertex_by_index( vtx );
      c[0] += mDerivs[3*i  ] * coords;
      c[1] += mDerivs[3*i+1] * coords;
      c[2] += mDerivs[3*i+2] * coords;
        // Remove data corresponding to fixed vertics from lists.
      if (vtx < pd.num_free_vertices()) {
        indices[w] = vtx;
        if (i != w) {
          mDerivs[3*w  ] = mDerivs[3*i  ];
          mDerivs[3*w+1] = mDerivs[3*i+1];
          mDerivs[3*w+2] = mDerivs[3*i+2];
        }
        ++w;
      }
    }
    indices.resize(w);
    mDerivs.resize(3*w);
    
    MsqMatrix<3,3> A( (MsqMatrix<3,1>*)c );
    

    MsqMatrix<1,3> tmp[3][3], gi;
    MsqMatrix<3,1> gj;
    MsqMatrix<3,3> W, dmdA, d2mdA2[6];
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    gradient.resize(w);
    Hessian.resize(w*(w+1)/2);
    size_t h = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
      gi = MsqMatrix<1,3>(&mDerivs[3*i]);
      gradient[i] = Vector3D( (dmdA * transpose(gi)).data() );
      tmp[0][0] = gi * d2mdA2[0];
      tmp[0][1] = gi * d2mdA2[1];
      tmp[0][2] = gi * d2mdA2[2];
      tmp[1][0] = gi * transpose(d2mdA2[1]);
      tmp[1][1] = gi * d2mdA2[3];
      tmp[1][2] = gi * d2mdA2[4];
      tmp[2][0] = gi * transpose(d2mdA2[2]);
      tmp[2][1] = gi * transpose(d2mdA2[4]);
      tmp[2][2] = gi * d2mdA2[5];
     
      Matrix3D& H = Hessian[h++];
      H[0][0] =           tmp[0][0] * transpose(gi);
      H[0][1] = H[1][0] = tmp[0][1] * transpose(gi);
      H[0][2] = H[2][0] = tmp[0][2] * transpose(gi);
      H[1][1] =           tmp[1][1] * transpose(gi);
      H[1][2] = H[2][1] = tmp[1][2] * transpose(gi);
      H[2][2] =           tmp[2][2] * transpose(gi);
      for (size_t j = i+1; j < indices.size(); ++j) {
        Matrix3D& H = Hessian[h++];
        gj = MsqMatrix<3,1>(&mDerivs[3*j]);
        H[0][0] = tmp[0][0] * gj;
        H[0][1] = tmp[0][1] * gj;
        H[0][2] = tmp[0][2] * gj;
        H[1][0] = tmp[1][0] * gj;
        H[1][1] = tmp[1][1] * gj;
        H[1][2] = tmp[1][2] * gj;
        H[2][0] = tmp[2][0] * gj;
        H[2][1] = tmp[2][1] * gj;
        H[2][2] = tmp[2][2] * gj;
      }
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
    for (size_t i = 0; i < gradient.size(); ++i)
      gradient[i] *= ck;
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
                                           msq_std::vector<Vector3D>& gradient,
                                           msq_std::vector<SymMatrix3D>& diagonal,
                                           MsqError& err )
{
    // make sure reinterpret_casts below are valid
  assert( sizeof(MsqMatrix<3,1>) == sizeof(Vector3D) );

  const unsigned s = ElemSampleQM::sample( handle );
  const size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  unsigned edim = TopologyInfo::dimension( elem.get_element_type() );
  const size_t* conn = elem.get_vertex_index_array();
  mapping_function_derivs( pd, handle, indices, mDerivs, err ); MSQ_ERRZERO(err);
  
  if (edim != 3) // use finite difference approximation for surface elements
    return QualityMetric::evaluate_with_Hessian_diagonal( pd, handle, value, indices, gradient, diagonal, err );
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for TMP metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    Vector3D c[3] = { Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0) };
    size_t w = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        // Convert from indices into element connectivity list to
        // indices into vertex array in patch data.
      size_t vtx = conn[indices[i]];
        // calculate Jacobian
      Vector3D coords = pd.vertex_by_index( vtx );
      c[0] += mDerivs[3*i  ] * coords;
      c[1] += mDerivs[3*i+1] * coords;
      c[2] += mDerivs[3*i+2] * coords;
        // Remove data corresponding to fixed vertics from lists.
      if (vtx < pd.num_free_vertices()) {
        indices[w] = vtx;
        if (i != w) {
          mDerivs[3*w  ] = mDerivs[3*i  ];
          mDerivs[3*w+1] = mDerivs[3*i+1];
          mDerivs[3*w+2] = mDerivs[3*i+2];
        }
        ++w;
      }
    }
    indices.resize(w);
    mDerivs.resize(3*w);
    
    MsqMatrix<3,3> A( (MsqMatrix<3,1>*)c );
    

    MsqMatrix<3,3> W, dmdA, d2mdA2[6];
    MsqMatrix<3,1> g;
    targetCalc->get_3D_target( pd, e, samplePts, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); MSQ_ERRZERO(err);
    gradient.resize(w);
    diagonal.resize(w);
    for (size_t i = 0; i < indices.size(); ++i) {
      g = MsqMatrix<3,1>(&mDerivs[3*i]);
      gradient[i] = Vector3D( (dmdA * g).data() );
      SymMatrix3D& H = diagonal[i];
      for (unsigned j = 0; j < 6; ++j)
        H[j] = transpose(g) * d2mdA2[j] * g; 
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
      gradient[i] *= ck;
      diagonal[i] *= ck;
    }
  }
  
  return rval;
}




} // namespace Mesquite
