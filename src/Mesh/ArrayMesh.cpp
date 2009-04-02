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


/** \file ArrayMesh.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ArrayMesh.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"

namespace MESQUITE_NS {

class IndexIterator : public EntityIterator
{
public:
  IndexIterator( size_t mStart, size_t mEnd )
    : mStart(mStart), mEnd(mEnd), mCurrent(mStart) {}
  virtual ~IndexIterator() {}
  virtual void restart() { mCurrent = mStart; }
  virtual Mesh::EntityHandle operator*() const
    { return (Mesh::EntityHandle)mCurrent; }
  virtual void operator++() { ++mCurrent; }
  virtual bool is_at_end() const { return mEnd - mCurrent <= 1; }
private:
  size_t mStart, mEnd, mCurrent;
};

ArrayMesh::ArrayMesh()
  : mDimension( 0 ),
    vertexCount( 0 ),
    coordArray( 0 ),
    fixedFlags( 0 ),
    vertexByteArray( 0 ),
    elementCount( 0 ),
    connArray( 0 ),
    connOffsets( 0 ),
    allocConnOffsets( 0 ),
    elementType( MIXED ),
    elementTypes( 0 ),
    nodesPerElement( 0 ),
    oneBasedArrays( false ),
    vertexAdjacencyList(0),
    vertexAdjacencyOffsets(0)
{}

ArrayMesh::ArrayMesh( int coords_per_vertex,
                      unsigned long num_vertices,
                      double* interleaved_vertex_coords,
                      const int* vertex_fixed_flags,
                      unsigned long num_elements,
                      EntityTopology element_type,
                      const unsigned long* element_connectivity_array,
                      bool one_based_conn_indices,
                      unsigned nodes_per_element ) 
  : mDimension( coords_per_vertex ),
    vertexCount( num_vertices ),
    coordArray( interleaved_vertex_coords ),
    fixedFlags( vertex_fixed_flags ),
    vertexByteArray( new unsigned char[num_vertices + one_based_conn_indices] ),
    elementCount( num_elements ),
    connArray( element_connectivity_array ),
    connOffsets( 0 ),
    allocConnOffsets( 0 ),
    elementType( element_type ),
    elementTypes( 0 ),
    nodesPerElement( nodes_per_element ),
    oneBasedArrays( one_based_conn_indices ),
    vertexAdjacencyList(0),
    vertexAdjacencyOffsets(0)
{
  if (oneBasedArrays) {
    coordArray -= mDimension;
    --fixedFlags;
  }
  
  if (nodesPerElement < 2)
    nodesPerElement = TopologyInfo::corners( element_type );
  memset( vertexByteArray, 0, num_vertices + one_based_conn_indices );
}

ArrayMesh::ArrayMesh( int coords_per_vertex,
                      unsigned long num_vertices,
                      double* interleaved_vertex_coords,
                      const int* vertex_fixed_flags,
                      unsigned long num_elements,
                      const EntityTopology* element_types,
                      const unsigned long* element_connectivity_array,
                      const unsigned long* element_connectivity_offsets,
                      bool one_based_conn_indices ) 
  : mDimension( coords_per_vertex ),
    vertexCount( num_vertices ),
    coordArray( interleaved_vertex_coords ),
    fixedFlags( vertex_fixed_flags ),
    vertexByteArray( new unsigned char[num_vertices + one_based_conn_indices] ),
    elementCount( num_elements ),
    connArray( element_connectivity_array ),
    connOffsets( element_connectivity_offsets ),
    allocConnOffsets( 0 ),
    elementType( MIXED ),
    elementTypes( element_types ),
    nodesPerElement( 0 ),
    oneBasedArrays( one_based_conn_indices ),
    vertexAdjacencyList(0),
    vertexAdjacencyOffsets(0)
{
  if (oneBasedArrays) {
    coordArray -= mDimension;
    --fixedFlags;
    if (element_connectivity_offsets)
      --connArray;
  }
  
  if (!element_connectivity_offsets) {
    connOffsets = allocConnOffsets = new unsigned long[num_elements+1];
    allocConnOffsets[0] = 0;
    for (unsigned long i = 1; i <= num_elements; ++i)
      allocConnOffsets[i] = allocConnOffsets[i-1] + TopologyInfo::corners( elementTypes[i-1] );
  }
  
  memset( vertexByteArray, 0, num_vertices + one_based_conn_indices );
}

void ArrayMesh::clear_mesh()
{
  delete [] vertexByteArray;
  delete [] vertexAdjacencyList;
  delete [] vertexAdjacencyOffsets;
  delete [] allocConnOffsets;
  mDimension = 0;
  vertexCount = 0;
  coordArray = 0;
  connOffsets = 0;
  allocConnOffsets = 0;
  fixedFlags = 0;
  vertexByteArray = 0;
  elementCount = 0;
  elementType = MIXED;
  elementTypes = 0;
  nodesPerElement = 0;
  oneBasedArrays = false;
  vertexAdjacencyList = 0;
  vertexAdjacencyOffsets = 0;
}

void ArrayMesh::set_mesh( int coords_per_vertex,
                          unsigned long num_vertices,
                          double* interleaved_vertex_coords,
                          const int* vertex_fixed_flags,
                          unsigned long num_elements,
                          EntityTopology element_type,
                          const unsigned long* element_connectivity_array,
                          bool one_based_conn_indices,
                          unsigned nodes_per_element ) 
{
  clear_mesh();
  mDimension = coords_per_vertex;
  vertexCount = num_vertices;
  coordArray = interleaved_vertex_coords;
  fixedFlags = vertex_fixed_flags;
  elementCount = num_elements;
  connArray = element_connectivity_array;
  elementType = element_type;
  oneBasedArrays = one_based_conn_indices;
  
  if (oneBasedArrays) {
    coordArray -= mDimension;
    --fixedFlags;
  }
  
  if (nodes_per_element < 2)
    nodesPerElement = TopologyInfo::corners( element_type );
  else
    nodesPerElement = nodes_per_element;
    
  vertexByteArray = new unsigned char[num_vertices + one_based_conn_indices];
  memset( vertexByteArray, 0, num_vertices + one_based_conn_indices );
}
  
  

ArrayMesh::~ArrayMesh()
{
  delete [] vertexByteArray;
  delete [] vertexAdjacencyList;
  delete [] vertexAdjacencyOffsets;
  delete [] allocConnOffsets;
}

inline const unsigned long* ArrayMesh::elem_verts( size_t e, int& n ) const
{
  if (connOffsets) {
    n = connOffsets[e+1] - connOffsets[e];
    return connArray + connOffsets[e];
  }
  else {
    n = nodesPerElement;
    return connArray + nodesPerElement*e;
  }
}

int ArrayMesh::get_geometric_dimension( MsqError& )
  { return mDimension; }

void ArrayMesh::get_all_elements( msq_std::vector<ElementHandle>& elements, MsqError& )
{
  elements.resize( elementCount );
  for (unsigned long i = 0; i < elementCount; ++i)
    elements[i] = (Mesh::ElementHandle)i;
}

void ArrayMesh::get_all_vertices( msq_std::vector<VertexHandle>& vertices, MsqError& )
{
  vertices.resize( vertexCount );
  for (unsigned long i = 0; i < vertexCount; ++i)
    vertices[i] = (Mesh::VertexHandle)(i+oneBasedArrays);
}

VertexIterator* ArrayMesh::vertex_iterator( MsqError& )
  { return new IndexIterator( oneBasedArrays, vertexCount + oneBasedArrays ); }

ElementIterator* ArrayMesh::element_iterator( MsqError& err )
  { return new IndexIterator( 0, elementCount ); }

void ArrayMesh::vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                         bool fixed_flag_array[],
                                         size_t num_vtx, 
                                         MsqError & )
{
  const size_t* indices = (const size_t*)vert_array;
  for (size_t i = 0; i < num_vtx; ++i)
    fixed_flag_array[i] = !!fixedFlags[indices[i]];
}

void ArrayMesh::vertices_get_slaved_flag( const VertexHandle*, 
                                          bool*,
                                          size_t , 
                                          MsqError &err )
{
  MSQ_SETERR(err)("ArrayMesh::vertices_get_slaved_flag not implemented.  "
                  "Cannot support Settings::SLAVE_FLAG", MsqError::NOT_IMPLEMENTED);
}

void ArrayMesh::vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx,
                                           MsqError &err )
{
  const size_t* indices = (const size_t*)vert_array;
  if (mDimension == 3) 
    for (size_t i = 0; i < num_vtx; ++i) 
      coordinates[i].set( coordArray+3*indices[i] );
  else if (mDimension == 2) 
    for (size_t i =0; i < num_vtx; ++i)
      coordinates[i].set( coordArray[2*indices[i]], coordArray[2*indices[i]+1], 0.0 );
  else
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
}

void ArrayMesh::vertex_set_coordinates( VertexHandle vert,
                                          const Vector3D& coordinates,
                                           MsqError &err )
{
  size_t i = (size_t)vert;
  if (mDimension == 3) 
    coordinates.get_coordinates(coordArray+3*i);
  else if (mDimension == 2) {
    coordArray[2*i] = coordinates[0];
    coordArray[2*i+1] = coordinates[1];
  }
  else
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
}

void ArrayMesh::vertex_set_byte( VertexHandle vertex, unsigned char byte, MsqError &)
{
  vertexByteArray[(size_t)vertex] = byte;
}

void ArrayMesh::vertices_set_byte( const VertexHandle *vert_array,
                                   const unsigned char *byte_array,
                                   size_t array_size, 
                                   MsqError &err )
{
  const size_t* indices = (const size_t*)vert_array;
  for (size_t i = 0; i < array_size; ++i)
    vertexByteArray[indices[i]] = byte_array[i];
}

void ArrayMesh::vertex_get_byte( VertexHandle vertex, unsigned char* byte, MsqError &)
{
  *byte = vertexByteArray[(size_t)vertex];
}

void ArrayMesh::vertices_get_byte( const VertexHandle *vert_array,
                                   unsigned char *byte_array,
                                   size_t array_size, 
                                   MsqError &err )
{
  const size_t* indices = (const size_t*)vert_array;
  for (size_t i = 0; i < array_size; ++i)
    byte_array[i] = vertexByteArray[indices[i]];
}

void ArrayMesh::vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         msq_std::vector<ElementHandle>& elements,
                         msq_std::vector<size_t>& offsets,
                         MsqError& )
{
  const size_t* indices = (const size_t*)vertex_array;
  if (!vertexAdjacencyList)
    build_vertex_adjacency_list();
  
  elements.clear();
  offsets.resize( num_vertex + 1 );
  for (size_t i = 0; i < num_vertex; ++i) {
    offsets[i] = elements.size();
    for (size_t j = vertexAdjacencyOffsets[indices[i]];
         j < vertexAdjacencyOffsets[indices[i]+1]; ++j)
      elements.push_back( (ElementHandle)vertexAdjacencyList[j] );
  }
  offsets[num_vertex] = elements.size();
}

void ArrayMesh::elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   msq_std::vector<VertexHandle>& vert_handles,
                                   msq_std::vector<size_t>& offsets, 
                                   MsqError &)
{
  const size_t* indices = (const size_t*)elem_handles;
  offsets.resize( num_elems + 1);
  vert_handles.clear();
  for (size_t i = 0; i < num_elems; ++i) {
    int count;
    const unsigned long* conn = elem_verts( indices[i], count );
    size_t prev_size = vert_handles.size();
    offsets[i] = prev_size;
    vert_handles.resize( prev_size + count );
    std::copy( conn, conn+count, (size_t*)(&vert_handles[prev_size]) );
  }
  offsets[num_elems] = vert_handles.size();
}

void ArrayMesh::elements_get_topologies( const ElementHandle *handles,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError& )
{
  const size_t* indices = (const size_t*)handles;
  if (elementType == MIXED) 
    for (size_t i = 0; i < num_elements; ++i)
      element_topologies[i] = elementTypes[indices[i]];
  else 
    for (size_t i = 0; i < num_elements; ++i)
      element_topologies[i] = elementType;
}

void ArrayMesh::release_entity_handles( const EntityHandle*, size_t, MsqError& err )
{
  MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED);
}

void ArrayMesh::release() {}

void ArrayMesh::build_vertex_adjacency_list()
{
  delete [] vertexAdjacencyList;
  delete [] vertexAdjacencyOffsets;
  vertexAdjacencyOffsets = new unsigned long[vertexCount+oneBasedArrays+1];
  
    // for each vertex, store the number of elements the previous
    // vertex occurs in.
  memset( vertexAdjacencyOffsets, 0, sizeof(unsigned long)*(vertexCount+oneBasedArrays+1) );
  for (size_t i = 0; i < elementCount; ++i) {
    int n;
    const unsigned long* conn = elem_verts( i, n );
    for (int j = 0; j < n; ++j)
      ++vertexAdjacencyOffsets[conn[j]+1];
  }
  
    // convert vertexAdjacencyOffsets from a shifted list of counts
    // to a list of offsts
  for (size_t i = 1; i <= vertexCount+oneBasedArrays; ++i)
    vertexAdjacencyOffsets[i] += vertexAdjacencyOffsets[i-1];
    
    // allocate space and populate with reverse connectivity
  vertexAdjacencyList = new unsigned long[vertexAdjacencyOffsets[vertexCount+oneBasedArrays]];
  for (size_t i = 0; i < elementCount; ++i) {
    int n;
    const unsigned long* conn = elem_verts( i, n );
    for (int j = 0; j < n; ++j)
      vertexAdjacencyList[vertexAdjacencyOffsets[conn[j]]++] = i;
  }
  
  for (size_t i = vertexCount+oneBasedArrays; i > 0; --i)
    vertexAdjacencyOffsets[i] = vertexAdjacencyOffsets[i-1];
  vertexAdjacencyOffsets[0] = 0; 
}

TagHandle ArrayMesh::tag_create( const msq_std::string&, TagType, unsigned, const void*, MsqError &err)
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); return (TagHandle)-1; }
void ArrayMesh::tag_delete( TagHandle, MsqError& err )
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); }
TagHandle ArrayMesh::tag_get( const msq_std::string&, MsqError& err )
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); return (TagHandle)-1; }
void ArrayMesh::tag_properties( TagHandle, msq_std::string&, TagType&, unsigned&, MsqError& err )
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); }
void ArrayMesh::tag_set_element_data( TagHandle, size_t, const ElementHandle*, const void*, MsqError& err )
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); }
void ArrayMesh::tag_set_vertex_data ( TagHandle, size_t, const VertexHandle*, const void*, MsqError& err )
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); }
void ArrayMesh::tag_get_element_data( TagHandle, size_t, const ElementHandle*, void*, MsqError& err )
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); }
void ArrayMesh::tag_get_vertex_data ( TagHandle, size_t, const VertexHandle*, void*, MsqError& err )
  { MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); }
  
} // namespace Mesquite
