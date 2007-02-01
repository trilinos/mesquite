/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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


/** \file VertexQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "VertexQM.hpp"
#include "PatchData.hpp"
#include "ElemSampleQM.hpp"

namespace Mesquite {

void VertexQM::get_evaluations( PatchData& pd, 
                                msq_std::vector<size_t>& handles,
                                bool /*free_vertices_only*/, 
                                MsqError& err )
{
  get_vertex_evaluations( pd, handles, true, err );
}

void VertexQM::get_vertex_evaluations( PatchData& pd, 
                                       msq_std::vector<size_t>& handles,
                                       bool free_vertices_only, 
                                       MsqError& err )
{
  size_t count = free_vertices_only ? pd.num_free_vertices() : pd.num_nodes();
  handles.resize(count);
  for (size_t i = 0; i < count; ++i)
    handles[i] = i;
}

void VertexQM::get_vertex_corner_handles( PatchData& pd, 
                                          size_t vtx_idx,
                                          msq_std::vector<size_t>& handles,
                                          MsqError& err )
{
  size_t len, *elems = pd.get_vertex_element_adjacencies( vtx_idx, len, err );
  MSQ_ERRRTN(err);
  
  handles.resize(len);
  for (size_t i = 0; i < len; ++i) {
    const MsqMeshEntity& elem = pd.element_by_index( elems[i] );
    const size_t* verts = elem.get_vertex_index_array();
    const size_t* ptr = std::find( verts, verts+elem.node_count(), vtx_idx );
    unsigned idx = ptr - verts;
    handles[i] = ElemSampleQM::handle( idx, elems[i] );
  }
}

} // namespace Mesquite