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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   MeshTSTT.hpp
  \brief  


  \author Thomas Leurent
  \date   2003-06-12
*/

#ifndef MESQUITE_MESH_TSTT_HPP
#define MESQUITE_MESH_TSTT_HPP

//kkc 040203 replace with more general TSTT interfaces #include "TSTT_LocalTSTTMesh.hh"
#include "TSTT_Mesh.hh"
#include "TSTT_CoreEntitySetQuery.hh"
#include "TSTT_AdvancedEntitySetQuery.hh"
#include "TSTT_Tag.hh"
#include "TSTT_ModifiableMesh.hh"

#include "MeshInterface.hpp"

//class TSTT::LocalTSTTMesh;

namespace Mesquite
{
  /*!  \class MeshTSTT

  \brief MeshTSTT is an Adapter of the TSTT interface to
         the Mesquite Mesh interface. 
    
         The Object Adapter pattern (see Gamma & al.) is used to
         use a TSTT implementation within the Mesquite framework.
         
         In essence, the MeshTSTT class takes a TSTT mesh handle
         in its constructor and plugs at compile-time a TSTT
         implementation to the Mesquite framework.
  */
  class MeshTSTT : public Mesquite::Mesh
  {
    typedef void* TagHandle;
  public:
//********* Functions that are NOT inherited ************
//kkc 040203    MeshTSTT(TSTT::LocalTSTTMesh& tstt_mesh, MsqError &err);
    MeshTSTT(TSTT::Mesh& tstt_mesh, MsqError &err);
    virtual ~MeshTSTT();
    
//********* Functions that ARE inherited ************
      // Returns whether this mesh lies in a 2D or 3D coordinate system.
    virtual int get_geometric_dimension(Mesquite::MsqError &/*err*/) const;
    
      // Returns the number of entities of the indicated type.
    virtual size_t get_total_vertex_count(MsqError &err) const;
    virtual size_t get_total_element_count(MsqError &err) const;
    
      // Fills array with handles to all vertices/elements
      // in the mesh.
    virtual void get_all_vertices(VertexHandle *vert_array,
                            size_t array_size, MsqError &err);
    virtual void get_all_elements(ElementHandle *elem_array,
                            size_t array_size, MsqError &err);
    
      // Returns a pointer to an iterator that iterates over the
      // set of all vertices in this mesh.  The calling code should
      // delete the returned iterator when it is finished with it.
      // If vertices are added or removed from the Mesh after obtaining
      // an iterator, the behavior of that iterator is undefined.
    virtual VertexIterator* vertex_iterator(MsqError &err);
    
      // Returns a pointer to an iterator that iterates over the
      // set of all top-level elements in this mesh.  The calling code should
      // delete the returned iterator when it is finished with it.
      // If elements are added or removed from the Mesh after obtaining
      // an iterator, the behavior of that iterator is undefined.
    virtual ElementIterator* element_iterator(MsqError &err);

//************ Vertex Properties ********************
      // Returns true or false, indicating whether the vertex
      // is allowed to be repositioned.  True indicates that the vertex
      // is fixed and cannot be moved.  Note that this is a read-only
      // property; this flag can't be modified by users of the
      // Mesquite::Mesh interface.
    virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err);

      // Returns true or false, indicating whether the vertex
      // is on the boundary.  Boundary nodes may be treated as
      // a special case by some algorithms or culling methods.
      // Note that this is a read-only
      // property; this flag can't be modified by users of the
      // Mesquite::Mesh interface.
    virtual void vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
                                 size_t num_vtx, MsqError &err);
    
      // Get/set location of a vertex
    virtual void vertices_get_coordinates(VertexHandle vert_array[],
                                  MsqVertex* const &coordinates,
                                  const size_t &num_vtx, MsqError &err);
    virtual void vertex_set_coordinates(VertexHandle vertex,
                                 const Vector3D &coordinates, MsqError &err);
    
      // Each vertex has a byte-sized flag that can be used to store
      // flags.  This byte's value is neither set nor used by the mesh
      // implementation.  It is intended to be used by Mesquite algorithms.
      // Until a vertex's byte has been explicitly set, its value is 0.
    virtual void vertex_set_byte (VertexHandle vertex,
                            unsigned char byte, MsqError &err);
    virtual void vertices_set_byte (VertexHandle *vert_array,
                              unsigned char *byte_array,
                              size_t array_size, MsqError &err);
    
      // Retrieve the byte value for the specified vertex or vertices.
      // The byte value is 0 if it has not yet been set via one of the
      // *_set_byte() functions.
    virtual void vertex_get_byte(VertexHandle vertex,
                                 unsigned char *byte, MsqError &err);
    virtual void vertices_get_byte(VertexHandle *vert_array,
                                   unsigned char *byte_array,
                                   size_t array_size, MsqError &err);
    
//**************** Vertex Topology *****************    
      // Gets the number of elements attached to this vertex.
      // Useful to determine how large the "elem_array" parameter
      // of the vertex_get_attached_elements() function must be.
    virtual size_t vertex_get_attached_element_count(VertexHandle vertex, MsqError &err) const;
    
      // Gets the elements attached to this vertex.
    virtual void vertex_get_attached_elements(VertexHandle vertex,
                                              ElementHandle* elem_array,
                                              size_t sizeof_elem_array,
                                              MsqError &err);
    
//*************** Element Topology *************
    
      // Gets the number of vertices in this element.
      // This data can also be found by querying the
      // element's topology and getting the number
      // of vertices per element for that topology type.
    virtual size_t element_get_attached_vertex_count(ElementHandle elem,
                                                  MsqError &err) const;
    
// Returns the vertices that are part of the topological definition of each
// element in the "elem_handles" array.  When this function is called, the
// following must be true:
//   a) "elem_handles" points at an array of "num_elems" element handles.
//   b) "vert_handles" points at an array of size "sizeof_vert_handles"
//   c) "csr_data" points at an array of size "sizeof_csr_data"
//   d) "csr_offsets" points at an array of size "num_elems+1"
//      
// When this function returns, adjacency information will be stored
// in csr format:
//    a) "vert_handles" stores handles to all vertices found in one
//       or more of the elements.  Each vertex appears only
//       once in "vert_handles", even if it is in multiple elements.
//    b) "sizeof_vert_handles" is set to the number of vertex
//       handles placed into "vert_handles".
//    c) "sizeof_csr_data" is set to the total number of vertex uses (for
//       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
//       the two triangles share some vertices).
//    c) "csr_offsets" is filled such that csr_offset[i] indicates the location
//       of entity i's first adjacency in "csr_data".  The number of vertices
//       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
//       reason, csr_offsets[num_elems] is set to the new value of
//       "sizeof_csr_data".
//    d) "csr_data" stores integer offsets which give the location of
//       each adjacency in the "vert_handles" array.
//
// As an example of how to use this data, you can get the handle of the first
// vertex in element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ]
//
// and the second vertex of element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ]
// 
    virtual void elements_get_attached_vertices(ElementHandle *elem_handles,
                                                size_t num_elems,
                                                VertexHandle *vert_handles,
                                                size_t &sizeof_vert_handles,
                                                size_t *csr_data,
                                                size_t &sizeof_csr_data,
                                                size_t *csr_offsets,
                                                MsqError &err);
    
      // Identifies the vertices attached to the elements by returning
      // each vertex's global index.  The vertex's global index indicates
      // where that vertex can be found in the array returned by
      // Mesh::get_all_vertices.
    void elements_get_attached_vertex_indices(Mesquite::Mesh::ElementHandle elems[],
                                              size_t num_elems,
                                              size_t index_array[],
                                              size_t array_size,
                                              size_t* offsets,
                                              MsqError &err);
    
      // Returns the topology of the given entity.
    virtual EntityTopology element_get_topology(ElementHandle entity_handle,
                                           MsqError &err) const;
      // Returns the topologies of the given entities.  The "entity_topologies"
      // array must be at least "num_elements" in size.
    virtual void elements_get_topologies(ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err);
    
//**************** Memory Management ****************
      // Tells the mesh that the client is finished with a given
      // entity handle.  
    virtual void release_entity_handles(EntityHandle *handle_array,
                                        size_t num_handles, MsqError &err);
    
      // Instead of deleting a Mesh when you think you are done,
      // call release().  In simple cases, the implementation could
      // just call the destructor.  More sophisticated implementations
      // may want to keep the Mesh object to live longer than Mesquite
      // is using it.
    virtual void release();

  private:
    //kkc 040203 replace in favor of several interfaces    mutable ::TSTT::LocalTSTTMesh tsttMesh;
    
    mutable ::TSTT::Mesh tsttMesh;
    mutable ::TSTT::CoreEntitySetQuery tsttCoreQuery;
    mutable ::TSTT::AdvancedEntitySetQuery tsttAdvQuery;
    mutable ::TSTT::Tag tsttTag;
    mutable ::TSTT::ModifiableMesh tsttModMesh;

    msq_stdc::size_t nbVertices;
    msq_stdc::size_t nbElements;

    mutable ::TSTT::EntityType elementType; //!< Whether the mesh elements are
                                   //!< regions (3D) or faces (2D)
    
    TagHandle fixedVertexTag;
    TagHandle boundaryVertexTag;
    TagHandle vertexByteTag;
    unsigned char* mZeros;
    
    mutable ::SIDL::array<EntityHandle> oneEntity; //!< array has one entry only.
    mutable ::SIDL::array<TagHandle> oneTagValue; //!< array has one entry only.
    mutable ::SIDL::array<int32_t> oneInt; //!< array has one entry only.
    mutable ::SIDL::array<TSTT::EntityTopology> oneTopo; //!< array has one entry only.
    ::SIDL::array<double> threeDoubles; //!< array has 3 entries only..

    //! This array usually has size 0 when passed to TSTT functions with size 0. 
    mutable ::SIDL::array<EntityHandle> cachedAdjEntArray;
    mutable ::SIDL::array<int32_t> adjCsrPt;
    mutable ::SIDL::array<int32_t> adjCsrDat;
    mutable EntityHandle cachedVertex;
 };
}

#endif // MESQUITE_MESH_TSTT_HPP
