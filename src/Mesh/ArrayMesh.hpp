/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ArrayMesh.hpp
 *  \brief Access mesh stored in arrays
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_ARRAY_MESH_HPP
#define MSQ_ARRAY_MESH_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

namespace Mesquite {

class ArrayMesh : public Mesh
{
  public:
  
    ArrayMesh( int coords_per_vertex,
               unsigned long num_vertices,
               double* interleaved_vertex_coords,
               const int* vertex_fixed_flags,
               unsigned long num_elements,
               EntityTopology element_type,
               const unsigned long* element_connectivity_array,
               bool one_based_conn_indices = false,
               unsigned nodes_per_element = 0 );
    
    ~ArrayMesh();
    
    virtual int get_geometric_dimension( MsqError& err );

    virtual void get_all_elements( msq_std::vector<ElementHandle>& elements,
                                   MsqError& err );
    virtual void get_all_vertices( msq_std::vector<VertexHandle>& vertices,
                                   MsqError& err );
    
    virtual VertexIterator* vertex_iterator(MsqError &err);
    virtual ElementIterator* element_iterator(MsqError &err);

    virtual void vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                          bool fixed_flag_array[],
                                          size_t num_vtx, 
                                          MsqError &err );
   virtual void vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx,
                                           MsqError &err );
    virtual void vertex_set_coordinates( VertexHandle vertex,
                                         const Vector3D &coordinates,
                                         MsqError &err );
    
    virtual void vertex_set_byte( VertexHandle vertex,
                                  unsigned char byte, 
                                  MsqError &err);
    virtual void vertices_set_byte( const VertexHandle *vert_array,
                                    const unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err );
    
    virtual void vertex_get_byte( const VertexHandle vertex,
                                  unsigned char *byte, 
                                  MsqError &err );
    virtual void vertices_get_byte( const VertexHandle *vertex,
                                    unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err );
    
    virtual void vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         msq_std::vector<ElementHandle>& elements,
                         msq_std::vector<size_t>& offsets,
                         MsqError& err );
    
    virtual void elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   msq_std::vector<VertexHandle>& vert_handles,
                                   msq_std::vector<size_t>& offsets, 
                                   MsqError &err);
    
    virtual void elements_get_topologies(const ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err);

    
    virtual TagHandle tag_create( const msq_std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err);
    
    virtual void tag_delete( TagHandle handle, MsqError& err );
    
    
    virtual TagHandle tag_get( const msq_std::string& name, 
                               MsqError& err );
     
    virtual void tag_properties( TagHandle handle,
                                 msq_std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err );
    
    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err );

    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err );
    
    
    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err );
    
    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err );

    
    
    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles, 
                                         MsqError &err);
    
    virtual void release();

  private:
    
    void build_vertex_adjacency_list();
    
    int mDimension;
    unsigned long vertexCount;
    double* coordArray;
    const int* fixedFlags;
    unsigned char* vertexByteArray;
    
    unsigned long elementCount;
    const unsigned long* connArray;
    EntityTopology elementType;
    unsigned nodesPerElement;
    bool oneBasedArrays;
    
    unsigned long* vertexAdjacencyList;
    unsigned long* vertexAdjacencyOffsets;
};



} // namespace Mesquite

#endif
