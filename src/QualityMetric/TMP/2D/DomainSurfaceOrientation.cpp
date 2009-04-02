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


/** \file DomainSurfaceOrientation.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "DomainSurfaceOrientation.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MappingFunction.hpp"
#include "WeightCalculator.hpp"
#include "ElementQM.hpp"

namespace MESQUITE_NS {
int DomainSurfaceOrientation::get_negate_flag( ) const { return 1; }

msq_std::string DomainSurfaceOrientation::get_name() const
  { return msq_std::string("DomainSurfaceOrientation"); }

void DomainSurfaceOrientation::get_evaluations( PatchData& pd,
                                      msq_std::vector<size_t>& handles,
                                      bool free,
                                      MsqError& err )
{
  handles.clear();
  msq_std::vector<size_t> elems;
  ElementQM::get_element_evaluations( pd, elems, free, err ); MSQ_ERRRTN(err);
  for (msq_std::vector<size_t>::iterator i = elems.begin(); i != elems.end(); ++i)
  {
    NodeSet samples = pd.get_samples( *i );
    EntityTopology type = pd.element_by_index( *i ).get_element_type();
    if (TopologyInfo::dimension(type) == 2) {
      for (unsigned j = 0; j < TopologyInfo::corners(type); ++j)
        if (samples.corner_node(j))
          handles.push_back( handle( Sample(0,j), *i ) );
      for (unsigned j = 0; j < TopologyInfo::edges(type); ++j)
        if (samples.mid_edge_node(j)) 
          handles.push_back( handle( Sample(1,j), *i ) );
      if (samples.mid_face_node(0))
        handles.push_back( handle( Sample(2,0), *i ) );
    }
  }
}

void DomainSurfaceOrientation::get_element_evaluations( PatchData& pd,
                                              size_t elem,
                                              msq_std::vector<size_t>& handles,
                                              MsqError& err )
{
  EntityTopology type = pd.element_by_index( elem ).get_element_type();
  handles.clear();
  if (TopologyInfo::dimension(type) == 3) 
    return;

  NodeSet samples = pd.get_samples( elem );
  for (unsigned j = 0; j < TopologyInfo::corners(type); ++j)
    if (samples.corner_node(j))
      handles.push_back( handle( Sample(0,j), elem ) );
  for (unsigned j = 0; j < TopologyInfo::edges(type); ++j)
    if (samples.mid_edge_node(j)) 
      handles.push_back( handle( Sample(1,j), elem ) );
  if (samples.mid_face_node(0))
    handles.push_back( handle( Sample(2,0), elem ) );
}

bool DomainSurfaceOrientation::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  size_t num_idx;
  return evaluate_with_indices( pd, handle, value, mIndices, num_idx, err );
}

bool DomainSurfaceOrientation::evaluate_with_indices( PatchData& pd,
                                            size_t handle,
                                            double& value,
                                            msq_std::vector<size_t>& indices,
                                            MsqError& err )
{
  size_t num_idx = 0;
  bool rval = evaluate_with_indices( pd, handle, value, mIndices, num_idx, err );
  indices.resize( num_idx );
  std::copy( mIndices, mIndices+num_idx, indices.begin() );
  return rval;
}

bool DomainSurfaceOrientation::evaluate_with_indices( PatchData& pd,
                                            size_t handle,
                                            double& value,
                                            size_t* indices,
                                            size_t& num_idx,
                                            MsqError& err )
{
  Sample s = ElemSampleQM::sample( handle );
  size_t e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  if (TopologyInfo::dimension( type ) != 2) {
    MSQ_SETERR(err)("Surface Orientation metric is valid only for surface elements", 
                    MsqError::UNSUPPORTED_ELEMENT );
    return false;
  }
  
  const MappingFunction2D* func = pd.get_mapping_function_2D( type );
  if (!func) {
    MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
    return false;
  }
  
  const NodeSet bits = pd.non_slave_node_set( e );
  
  Vector3D n;
  pd.get_domain_normal_at_sample( e, s, n, err ); MSQ_ERRZERO(err);
  func->derivatives( s, bits, indices, mDerivs, num_idx, err ); MSQ_ERRZERO(err);
  func->convert_connectivity_indices( elem.node_count(), indices, num_idx, err ); MSQ_ERRZERO(err);
  
    // Convert from indices into element connectivity list to
    // indices into vertex array in patch data.
  const size_t* conn = elem.get_vertex_index_array();
  for (size_t* i = mIndices; i != mIndices+num_idx; ++i)
    *i = conn[*i];
  
  MsqVector<2>* d = mDerivs;
  Vector3D c[2] = { Vector3D(0,0,0), Vector3D(0,0,0) };
  for (size_t i = 0; i < num_idx; ++i, ++d) {
    Vector3D coords = pd.vertex_by_index( mIndices[i] );
    c[0] += (*d)[0] * coords;
    c[1] += (*d)[1] * coords;
  }
  
  n.normalize();
  const double m[2] = { c[0] % n, c[1] % n };
  value = m[0]*m[0] + m[1]*m[1];
  
    // apply target weight to value
  if (weightCalc) {
    const double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }
  return true;
}

} // namespace Mesquite
