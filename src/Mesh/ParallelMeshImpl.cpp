/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu
   
  ***************************************************************** */
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 26-March-08 at 10:26:21
//  LAST-MOD: 15-Nov-04 by kraftche@cae.wisc.edu
//
/*! \file ParallelMeshImpl.cpp

\brief This files contains a parallel mesh implementation that can be used
to run mesquite by default.

    \author Jason Kraftcheck
    \author Martin Isenburg
    \date 2008-3-26
 */

#include "ParallelMeshImpl.hpp"
#include "MeshImplData.hpp"
#include "MeshImplTags.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"

namespace Mesquite
{

ParallelMeshImpl::ParallelMeshImpl()
{
  helper = 0;
}

ParallelMeshImpl::ParallelMeshImpl(int num_nodes, int num_elem, 
                   Mesquite::EntityTopology entity_topology, 
                   double **coords, int **connectivity)
{
  int i,j;
  
  helper = 0;
  numCoords = 3;
  
  /*  
  std::cout << "in new mesh impl constructor "<< num_nodes << " " << num_elem <<  std::endl;

  for (i=0;i<num_nodes;i++) {
    std::cout << "  Vtx " << i << ": " << coords[i][0] << " " << coords[i][1] << " " << coords[i][2] << std::endl;
  }

  if (entity_topology == Mesquite::TRIANGLE) {
    for (i=0;i<num_elem;i++) {
      std::cout << "  Conn " << i << ": " << connectivity[i][0] << " " << connectivity[i][1] << " " << connectivity[i][2] << std::endl;
    }
  } else {
    for (i=0;i<num_elem;i++) {
      std::cout << "  Conn " << i << ": " << connectivity[i][0] << " " << connectivity[i][1] << " " << connectivity[i][2] <<" " << connectivity[i][3]<< std::endl;
    }
  }
  */

  Mesquite::MsqError err;
  myMesh->allocate_vertices( num_nodes, err ); MSQ_ERRRTN(err);
  myMesh->allocate_elements( num_elem, err ); MSQ_ERRRTN(err);

  // Fill in the data
  for (i = 0; i < num_nodes; ++i) {
     myMesh->reset_vertex( i, Vector3D(coords[i][0],coords[i][1],coords[i][2]), 
                           false, err) ;
  }

  int verts_per_elem=0;
  EntityTopology elem_type;
  
  if (entity_topology == Mesquite::TRIANGLE) {
    verts_per_elem = 3;
    elem_type = Mesquite::TRIANGLE;
  } else if (entity_topology == Mesquite::QUADRILATERAL) {
    verts_per_elem = 4;
    elem_type = Mesquite::QUADRILATERAL;
  } else if (entity_topology == Mesquite::TETRAHEDRON) {
    verts_per_elem = 4;
    elem_type = Mesquite::TETRAHEDRON;
  } else {
    std::cout << "error... only supporting triangles, quadrangles, and tets at this time "<< std::endl;
  }

  msq_std::vector<int> conn;
  
  if (conn.size() < (unsigned)(num_elem*verts_per_elem))
    conn.resize( num_elem*verts_per_elem );

  msq_std::vector<size_t> connect(verts_per_elem);
  for (i = 0; i < num_elem; i++) {
    for (j = 0; j < verts_per_elem; j++) {
      connect[j] = connectivity[i][j];
    }
    myMesh->reset_element( i, connect, elem_type, err);
  }
  
/** Get TAGS needed for parallel operations or create them if they don't
   already exist */

  const char GLOBAL_ID_NAME[] = "GLOBAL_ID";
  const char PROCESSOR_ID_NAME[] = "PROCESSOR_ID";
    
  gid_tag = tag_get( GLOBAL_ID_NAME, err );
    // if tag doesn't exist, create it
  if (!gid_tag) {
    err.clear();

    int default_gid[0];
    default_gid[0] = -1;

    gid_tag = tag_create( GLOBAL_ID_NAME, INT, 1, default_gid, err );
      // the 'INT' is the type of data to store
      // the '1' is for one value per vertex
      // NULL for no default value, if you want them all
      // initialized to something, pass in a pointer to an int
      // with the value.
  } 

  pid_tag = tag_get( PROCESSOR_ID_NAME, err );
    // if tag doesn't exist, create it
  if (!pid_tag) {
    err.clear();

    int default_pid[0];
    default_pid[0] = -1;
    pid_tag = tag_create( PROCESSOR_ID_NAME, INT, 1, default_pid, err );
  }
}

//**************** Parallel Methods ******************************

void ParallelMeshImpl::vertices_get_global_id(const VertexHandle vert_array[],
					      int gid[],
					      size_t num_vtx,
					      MsqError& err)
{
  tag_get_vertex_data( gid_tag, num_vtx, vert_array, gid, err );
  MSQ_CHKERR(err);
}

void ParallelMeshImpl::vertices_set_global_id(const VertexHandle vert_array[],
					      int gid[],
					      size_t num_vtx,
					      MsqError& err)
{
  tag_set_vertex_data( gid_tag, num_vtx, vert_array, gid, err );
  MSQ_CHKERR(err);
}

     
void ParallelMeshImpl::vertices_get_processor_id(const VertexHandle vert_array[],
						 int pid[],
						 size_t num_vtx,
						 MsqError& err)
{
  tag_get_vertex_data( pid_tag, num_vtx, vert_array, pid, err );
  MSQ_CHKERR(err);
}

void ParallelMeshImpl::vertices_set_processor_id(const VertexHandle vert_array[],
						 int pid[],
						 size_t num_vtx,
						 MsqError& err)
{
  tag_set_vertex_data( pid_tag, num_vtx, vert_array, pid, err );
  MSQ_CHKERR(err); 
}

void ParallelMeshImpl::set_parallel_helper(ParallelHelper* helper) {
  this->helper = helper;
}

ParallelHelper* ParallelMeshImpl::get_parallel_helper() {
  return helper;
}

} // namespace Mesquite

