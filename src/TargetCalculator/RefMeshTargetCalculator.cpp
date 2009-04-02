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


/** \file RefMeshTargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "ReferenceMesh.hpp"
#include "ElemSampleQM.hpp"

namespace MESQUITE_NS {

RefMeshTargetCalculator::~RefMeshTargetCalculator() {}

bool RefMeshTargetCalculator::get_3D_target( PatchData& pd, 
                                             size_t element,
                                             Sample sample,
                                             MsqMatrix<3,3>& W_out,
                                             MsqError& err )
{
  NodeSet ho_bits = get_vertex_coords( pd, element, err ); MSQ_ERRZERO(err);
  MsqMeshEntity& elem = pd.element_by_index( element );
  EntityTopology type = elem.get_element_type();
  const MappingFunction3D* func = pd.get_mapping_function_3D( type );
  jacobianCalc.get_Jacobian_3D( func, ho_bits, sample, &tmpCoords[0], 
                                elem.node_count(), W_out, err ); MSQ_ERRZERO(err);
  return true;
}

bool RefMeshTargetCalculator::get_2D_target( PatchData& pd, 
                                             size_t element,
                                             Sample sample,
                                             MsqMatrix<3,2>& W_out,
                                             MsqError& err )
{
  NodeSet ho_bits = get_vertex_coords( pd, element, err ); MSQ_ERRFALSE(err);
  MsqMeshEntity& elem = pd.element_by_index( element );
  EntityTopology type = elem.get_element_type();
  const MappingFunction2D* func = pd.get_mapping_function_2D( type );
  jacobianCalc.get_Jacobian_2D( func, ho_bits, sample, &tmpCoords[0], 
                                elem.node_count(), W_out, err ); MSQ_ERRZERO(err);
  return true;
}

NodeSet RefMeshTargetCalculator::get_vertex_coords( PatchData& pd, 
                                                     size_t elem_idx,  
                                                     MsqError& err )
{
  NodeSet bits;
  MsqMeshEntity& elem = pd.element_by_index( elem_idx );
  const EntityTopology type = elem.get_element_type();
  
  const unsigned n = elem.node_count();
  tmpCoords.resize( n );
  tmpHandles.resize( n );
  
  const msq_std::size_t* vtx_idx = elem.get_vertex_index_array();
  const Mesh::VertexHandle* vtx_hdl = pd.get_vertex_handles_array();
  
  for (unsigned i = 0; i < n; ++i)
    tmpHandles[i] = vtx_hdl[vtx_idx[i]];
  
  refMesh->get_reference_vertex_coordinates( &tmpHandles[0], n, &tmpCoords[0], err );
  if (MSQ_CHKERR(err)) return NodeSet();

  bool midedge, midface, midvol;
  TopologyInfo::higher_order( type, n, midedge, midface, midvol, err );
  if (MSQ_CHKERR(err)) return NodeSet();
  
  if (midedge) 
    bits.set_all_mid_edge_nodes(type);
  if (midface)
    bits.set_all_mid_face_nodes(type);
  if (TopologyInfo::dimension(type) == 3 && midvol)
    bits.set_mid_region_node();
 
  return bits;
}

} // namespace Mesquite
