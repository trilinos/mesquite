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


/** \file RefMeshCornerTarget.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "RefMeshCornerTarget.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "ReferenceMesh.hpp"
#include "MsqMatrix.hpp"
#include "ElemSampleQM.hpp"

namespace Mesquite {

RefMeshCornerTarget::~RefMeshCornerTarget() {}

bool RefMeshCornerTarget::get_3D_target( PatchData& pd, 
                                         size_t element,
                                         const SamplePoints* pts,
                                         Sample sample,
                                         MsqMatrix<3,3>& W_out,
                                         MsqError& err )
{
  EntityTopology type = pd.element_by_index(element).get_element_type();
  unsigned corner = sample.number;
  if (sample.dimension != 0) {
    MSQ_SETERR(err)("RefMeshCornerTarget cannot generate targets at sample "
                    "locations other than corners", MsqError::INVALID_STATE);
    return false;
  }
  
  Vector3D cols[3];
  unsigned num_adj;
  const unsigned* adj = TopologyInfo::adjacent_vertices( type, corner, num_adj );
  const Mesh::VertexHandle* handles = pd.get_vertex_handles_array();
  const size_t* conn = pd.element_by_index(element).get_vertex_index_array();
  
  vertHandles[0] = handles[conn[corner]];
  vertHandles[1] = handles[conn[adj[0]]];
  vertHandles[2] = handles[conn[adj[1]]];
  switch (num_adj) {
    case 3:
      vertHandles[3] = handles[conn[adj[2]]];
      refMesh->get_reference_vertex_coordinates( vertHandles, 4, refCoords, err );
      MSQ_ERRZERO(err);
      cols[0] = refCoords[1] - refCoords[0];
      cols[1] = refCoords[2] - refCoords[0];
      cols[2] = refCoords[3] - refCoords[0];
      break;
    case 2:
      refMesh->get_reference_vertex_coordinates( vertHandles, 3, refCoords, err );
      MSQ_ERRZERO(err);
      cols[0] = refCoords[1] - refCoords[0];
      cols[1] = refCoords[2] - refCoords[0];
      if (pd.domain_set()) {
        pd.get_domain_normal_at_corner( element, corner, cols[2], err );
        MSQ_ERRZERO(err);
      }
      else {
        cols[2] = cols[0] * cols[1];
      }
      cols[2].normalize();
      break;
    default:
      // pyramid apex or polyheron?
      MSQ_SETERR(err)("Cannot evaluate RefMeshCornerTarget at non-corner vertex.",
                      MsqError::UNSUPPORTED_ELEMENT);
      return false;
  }
  
  W_out(0,0) = cols[0][0];
  W_out(0,1) = cols[1][0];
  W_out(0,2) = cols[2][0];
  W_out(1,0) = cols[0][1];
  W_out(1,1) = cols[1][1];
  W_out(1,2) = cols[2][1];
  W_out(2,0) = cols[0][2];
  W_out(2,1) = cols[1][2];
  W_out(2,2) = cols[2][2];
  
  return true;
}

bool RefMeshCornerTarget::get_2D_target( PatchData& , 
                                         size_t ,
                                         const SamplePoints* ,
                                         Sample ,
                                         MsqMatrix<3,2>& ,
                                         MsqError& err )
{
  MSQ_SETERR(err)("RefMeshCornerTarget cannot be used with jacobian-based "
                  "target metrics", MsqError::INVALID_STATE);
  return false;
}

bool RefMeshCornerTarget::surface_targets_are_3D() const
  { return true; }

} // namespace Mesquite
