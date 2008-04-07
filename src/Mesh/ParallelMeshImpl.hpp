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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov  ,
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */
/*!
  \file   ParallelMeshImpl.hpp
  \brief  

  \author Darryl Melander
  \author Thomas Leurent
  \author Jason Kraftcheck
  \author Martin Isenburg
  \date   2003-04-17
*/

#ifndef PARALLEL_MESQUITE_MESH_IMPL_HPP
#define PARALLEL_MESQUITE_MESH_IMPL_HPP

#include "ParallelMeshInterface.hpp"
#include "MeshImpl.hpp"

namespace Mesquite
{
  /*!  \class ParallelMeshImpl
  \brief ParallelMeshImpl is a Mesquite implementation of the ParallelMesh
  interface. It inherits all of the implementation from MeshImpl and only
  implements any additional functionality.
  */

  class MESQUITE_EXPORT ParallelMeshImpl : public MeshImpl, virtual public ParallelMesh
  {
  public:
//********* Functions that are NOT inherited ************

    ParallelMeshImpl();

    ParallelMeshImpl(int num_vertex, int num_elem, 
             EntityTopology entity_topology, 
	     double **coords, int **conn);

//**************** Parallel Methods ******************************

    /*! Get/Set global ids for given vertices.
     */
    virtual void vertices_get_global_id(const VertexHandle vert_array[],
                                        int gid[],
                                        size_t num_vtx,
                                        MsqError& err);
     
    virtual void vertices_set_global_id(const VertexHandle vert_array[],
                                       int gid[],
                                       size_t num_vtx,
                                       MsqError& err);
     
    /*! Get/Set processor ids for given vertices.
     */
    virtual void vertices_get_processor_id(const VertexHandle vert_array[],
                                           int pid[],
                                           size_t num_vtx,
                                           MsqError& err);
     
    virtual void vertices_set_processor_id(const VertexHandle vert_array[],
                                           int pid[],
                                           size_t num_vtx,
                                           MsqError& err);
     
    /*! Get/Set parallel helper
     */
    virtual void set_parallel_helper(ParallelHelper* helper);
    virtual ParallelHelper* get_parallel_helper();

  private:
    void* gid_tag;
    void* pid_tag;
    ParallelHelper* helper;
  };
}

#endif
