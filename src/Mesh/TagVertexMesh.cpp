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


/** \file TagVertexMesh.cpp
 *  \brief Implementation of Mesquite::TagVertexMesh class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TagVertexMesh.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"

namespace Mesquite {

void TagVertexMesh::initialize( Mesh* mesh, msq_std::string name, bool init, MsqError& err )
{
  realMesh = mesh;
  tagName = name;

  tagHandle = get_mesh()->tag_get( tagName, err );
    // If tag isn't defined yet, we're done for now.
  if (err.error_code() == MsqError::TAG_NOT_FOUND) {
    err.clear();
    return;
  } 
  else if (MSQ_CHKERR(err))
    return;
  
    // If tag is already defined, make sure it is the correct type.
  msq_std::string t_name;
  Mesh::TagType type;
  unsigned length;
  tag_properties( tagHandle, t_name, type, length, err ); 
  MSQ_ERRRTN(err);
  if (!(type == Mesh::DOUBLE && length == 3) &&
      !(type == Mesh::BYTE && length == 3*sizeof(double)))
    MSQ_SETERR(err)(MsqError::TAG_ALREADY_EXISTS,
                    "Tag \"%s\" has invalid type or size.",
                    tagName.c_str());
 
    // If tag is already defined and init was true, reset tag
    // values.
  haveTagHandle = true;
  if (init) {
    copy_all_coordinates( err ); 
    MSQ_ERRRTN(err);
  }
}

void TagVertexMesh::copy_all_coordinates( MsqError& err )
{
  if (!haveTagHandle) {
    MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
    return;
  }

  msq_std::vector<Mesh::VertexHandle> handles;
  get_all_vertices( handles, err ); 
  MSQ_ERRRTN(err);
  
  msq_std::vector<MsqVertex> coords(handles.size());
  get_mesh()->vertices_get_coordinates( &handles[0], &coords[0], handles.size(), err );
  MSQ_ERRRTN(err);
  
  msq_std::vector<double> data( 3*handles.size() );
  msq_std::vector<double>::iterator j = data.begin();
  msq_std::vector<MsqVertex>::const_iterator i = coords.begin();
  while (i != coords.end()) {
    i->get_coordinates( &*j );
    ++i;
    j += 3;
  }
  
  tag_set_vertex_data( tagHandle, handles.size(), &handles[0], &data[0], err );
  MSQ_ERRRTN(err);
}

void TagVertexMesh::check_remove_tag( MsqError& err )
{
  if (cleanUpTag && haveTagHandle) {
    tag_delete( tagHandle, err );
    MSQ_ERRRTN(err);
  }
  haveTagHandle = false;
}
    

TagVertexMesh::TagVertexMesh( MsqError& err,
                              Mesh* real_mesh,
                              bool init,
                              bool clean_up,
                              msq_std::string name)
  : tagHandle(0), haveTagHandle(false), cleanUpTag(clean_up)
{
  if (name.size() == 0)
    name = "MsqAltCoords";
  initialize( real_mesh, name, init, err ); MSQ_CHKERR(err);
}

TagVertexMesh::~TagVertexMesh()
{
  MsqError err;
  check_remove_tag( err );
}

void TagVertexMesh::set_mesh( Mesh* mesh, bool init, MsqError& err )
{
  check_remove_tag( err ); MSQ_ERRRTN(err);
  initialize( mesh, tagName, init, err ); MSQ_ERRRTN(err);
}

void TagVertexMesh::set_tag_name( msq_std::string name, bool init, MsqError& err )
{
  check_remove_tag( err ); MSQ_ERRRTN(err);
  initialize( realMesh, name, init, err ); MSQ_ERRRTN(err);
}

void TagVertexMesh::clear( MsqError& err )
{
  if (haveTagHandle) {
    copy_all_coordinates( err );
    MSQ_CHKERR(err);
  }
}


void TagVertexMesh::vertices_get_coordinates( const VertexHandle vert_array[],
                                              MsqVertex* coordinates,
                                              size_t num_vtx,
                                              MsqError &err )
{
  if (!haveTagHandle) {
    get_mesh()->vertices_get_coordinates( vert_array, coordinates, num_vtx, err );
    MSQ_ERRRTN(err);
  }
  else {
    msq_std::vector<double> coords( num_vtx * 3 );
    get_mesh()->tag_get_vertex_data( tagHandle, num_vtx, vert_array, &coords[0], err );
    MSQ_ERRRTN(err);
    MsqVertex* coordinates_end = coordinates + num_vtx;
    msq_std::vector<double>::const_iterator i = coords.begin();
    while (coordinates != coordinates_end) {
      coordinates->set( &*i );
      i += 3;
      ++coordinates;
    }
  }
}
    
void TagVertexMesh::vertex_set_coordinates( VertexHandle vertex,
                                            const Vector3D &coordinates,
                                            MsqError &err )
{
  if (!haveTagHandle)
  {
    tagHandle = get_mesh()->tag_create( tagName, Mesh::DOUBLE, 3, 0, err ); 
    MSQ_ERRRTN(err);
    haveTagHandle = true;
    copy_all_coordinates( err );
    MSQ_ERRRTN(err);  
  }
  
  get_mesh()->tag_set_vertex_data( tagHandle, 1, &vertex, coordinates.to_array(), err );
  MSQ_ERRRTN(err);
}



//************ Operations on entire mesh ****************

int TagVertexMesh::get_geometric_dimension( MsqError& err )
  { return get_mesh()->get_geometric_dimension( err ); }

void TagVertexMesh::get_all_elements( msq_std::vector<ElementHandle>& elements, MsqError& err )
  { get_mesh()->get_all_elements( elements, err ); }

void TagVertexMesh::get_all_vertices( msq_std::vector<VertexHandle>& vertices, MsqError& err )
  { get_mesh()->get_all_vertices( vertices, err ); }


//************ Vertex Properties ********************

void TagVertexMesh::vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                             bool fixed_flag_array[],
                                             size_t num_vtx, 
                                             MsqError &err )
  { get_mesh()->vertices_get_fixed_flag( vert_array, fixed_flag_array, num_vtx, err ); }

void TagVertexMesh::vertex_set_byte( VertexHandle vertex,
                                     unsigned char byte, 
                                     MsqError &err)
  { get_mesh()->vertex_set_byte( vertex, byte, err ); }

void TagVertexMesh::vertices_set_byte( const VertexHandle *vert_array,
                                       const unsigned char *byte_array,
                                       size_t array_size, 
                                       MsqError &err )
  { get_mesh()->vertices_set_byte( vert_array, byte_array, array_size, err ); }

void TagVertexMesh::vertex_get_byte( const VertexHandle vertex,
                                     unsigned char *byte, 
                                     MsqError &err )
  { get_mesh()->vertex_get_byte( vertex, byte, err ); }

void TagVertexMesh::vertices_get_byte( const VertexHandle *vertex,
                                       unsigned char *byte_array,
                                       size_t array_size, 
                                       MsqError &err )
  { get_mesh()->vertices_get_byte( vertex, byte_array, array_size, err ); }
    
//**************** Vertex Topology *****************    

void TagVertexMesh::vertices_get_attached_elements( 
                                const VertexHandle* vertex_array,
                                size_t num_vertex,
                                msq_std::vector<ElementHandle>& elements,
                                msq_std::vector<size_t>& offsets,
                                MsqError& err )
  { get_mesh()->vertices_get_attached_elements( vertex_array, num_vertex, elements, offsets, err ); }
    
//*************** Element Topology *************

void TagVertexMesh::elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   msq_std::vector<VertexHandle>& vert_handles,
                                   msq_std::vector<size_t>& offsets, 
                                   MsqError &err)
  { get_mesh()->elements_get_attached_vertices( elem_handles, num_elems, vert_handles, offsets, err ); }

void TagVertexMesh::elements_get_topologies(
                                    const ElementHandle *element_handle_array,
                                    EntityTopology *element_topologies,
                                    size_t num_elements, MsqError &err)
  { get_mesh()->elements_get_topologies( element_handle_array, element_topologies, num_elements, err ); }

//***************  Tags  ***********

TagHandle TagVertexMesh::tag_create( const msq_std::string& tag_name,
                                     TagType type, unsigned length,
                                     const void* default_value,
                                     MsqError &err)
{
    // Don't allow access to internal tag for vertex coordinates.
    // This prevents accidental layering of multiple instances of
    // TagVertexMesh with the same tag name.
  if (tag_name == tagName) {
    MSQ_SETERR(err)("Attempt to access internal tag data using tag interface.",
                    MsqError::TAG_ALREADY_EXISTS);
    return (TagHandle)0;
  }
  
  return get_mesh()->tag_create( tag_name, type, length, default_value, err );
}

void TagVertexMesh::tag_delete( TagHandle handle, MsqError& err )
  { get_mesh()->tag_delete( handle, err ); }

TagHandle TagVertexMesh::tag_get( const msq_std::string& name, MsqError& err )
{
    // Don't allow access to internal tag for vertex coordinates.
    // This prevents accidental layering of multiple instances of
    // TagVertexMesh with the same tag name.
  if (name == tagName) {
    MSQ_SETERR(err)("Attempt to access internal tag data using tag interface.",
                    MsqError::INVALID_ARG);
    return (TagHandle)0;
  }
  
  return get_mesh()->tag_get( name, err );
}


void TagVertexMesh::tag_properties( TagHandle handle,
                                    msq_std::string& name_out,
                                    TagType& type_out,
                                    unsigned& length_out,
                                    MsqError& err )
  { get_mesh()->tag_properties( handle, name_out, type_out, length_out, err ); }

void TagVertexMesh::tag_set_element_data( TagHandle handle,
                                          size_t num_elems,
                                          const ElementHandle* elem_array,
                                          const void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_set_element_data( handle, num_elems, elem_array, tag_data, err ); }

void TagVertexMesh::tag_set_vertex_data ( TagHandle handle,
                                          size_t num_elems,
                                          const VertexHandle* node_array,
                                          const void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_set_vertex_data( handle, num_elems, node_array, tag_data, err ); }

void TagVertexMesh::tag_get_element_data( TagHandle handle,
                                          size_t num_elems,
                                          const ElementHandle* elem_array,
                                          void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_get_element_data( handle, num_elems, elem_array, tag_data, err ); }

void TagVertexMesh::tag_get_vertex_data ( TagHandle handle,
                                          size_t num_elems,
                                          const VertexHandle* node_array,
                                          void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_get_vertex_data( handle, num_elems, node_array, tag_data, err ); }
    
//**************** Memory Management ****************

void TagVertexMesh::release_entity_handles( const EntityHandle *handle_array,
                                            size_t num_handles, 
                                            MsqError &err)
  { get_mesh()->release_entity_handles( handle_array, num_handles, err ); }

void TagVertexMesh::release()
{
  MsqError err;
  clear(err);
  check_remove_tag(err);
  realMesh->release();
  realMesh = 0;
  haveTagHandle = false;
}


} // namespace Mesquite
