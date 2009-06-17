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

namespace MESQUITE_NS {


msq_std::string TagVertexMesh::get_name() const
{
  std::string result("TagVertexMesh(\"");
  result += tagName;
  result += "\")";
  return result;
}

void TagVertexMesh::initialize( Mesh* mesh, msq_std::string name, MsqError& err )
{
  MeshDecorator::set_mesh( mesh );
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

  haveTagHandle = true;
}


double TagVertexMesh::loop_over_mesh( Mesh* mesh, 
                                      MeshDomain* , 
                                      const Settings* ,
                                      MsqError& err )
{
  if (mesh != get_mesh()) {
    MSQ_SETERR(err)("InstructionQueue and TagVertexMesh have different "
                    "Mesquite::Mesh instances.  Cannot initialize TagVertexMesh",
                    MsqError::INVALID_MESH);
    return 0.0;
  }
  
  copy_all_coordinates( err ); MSQ_ERRZERO(err);
  return 0.0;
}

void TagVertexMesh::copy_all_coordinates( MsqError& err )
{
  if (!haveTagHandle) {
    tagHandle = get_mesh()->tag_create( tagName, Mesh::DOUBLE, 3, 0, err ); 
    MSQ_ERRRTN(err);
    haveTagHandle = true;
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
                              bool clean_up,
                              msq_std::string name)
  : tagHandle(0), haveTagHandle(false), cleanUpTag(clean_up)
{
  if (name.size() == 0)
    name = "MsqAltCoords";
  initialize( real_mesh, name, err ); MSQ_CHKERR(err);
}

TagVertexMesh::~TagVertexMesh()
{
  MsqError err;
  check_remove_tag( err );
}

void TagVertexMesh::set_mesh( Mesh* mesh, MsqError& err )
{
  check_remove_tag( err ); MSQ_ERRRTN(err);
  initialize( mesh, tagName, err ); MSQ_ERRRTN(err);
}

void TagVertexMesh::set_tag_name( msq_std::string name, MsqError& err )
{
  check_remove_tag( err ); MSQ_ERRRTN(err);
  initialize( get_mesh(), name, err ); MSQ_ERRRTN(err);
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
    copy_all_coordinates( err );
    MSQ_ERRRTN(err);  
  }
  
  get_mesh()->tag_set_vertex_data( tagHandle, 1, &vertex, coordinates.to_array(), err );
  MSQ_ERRRTN(err);
}





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
    
void TagVertexMesh::release()
{
  MsqError err;
  clear(err);
  check_remove_tag(err);
  haveTagHandle = false;
  MeshDecorator::release();
}


} // namespace Mesquite
