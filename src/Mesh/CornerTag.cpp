/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
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

    kraftche@cae.wisc.edu    

  ***************************************************************** */

#include "CornerTag.hpp"
#include "MeshSet.hpp"

namespace Mesquite {

size_t CornerTagHandles::size( Mesh::TagType type )
{
  switch (type) {
    case Mesh::BYTE:   return 1;
    case Mesh::BOOL:   return sizeof(bool);
    case Mesh::INT:    return sizeof(int);
    case Mesh::DOUBLE: return sizeof(double);
    case Mesh::HANDLE: return sizeof(Mesh::EntityHandle);
    default:           return 0;
  }
}

Mesh* CornerTagHandles::get_current_mesh( PatchData* pd )
{
  return pd->get_mesh_set() ? pd->get_mesh_set()->get_current_mesh() : 0;
}

TagHandle CornerTagHandles::get_handle( Mesh* mesh, unsigned corners, MsqError& err )
{
    // Resize vector as necessary
  if (corners >= cornerHandles.size())
  {
    cornerHandles.resize( corners+1, 0 );
  }
  
    // Get handle if we don't already have it
  if (!cornerHandles[corners])
  {
      // Construct tag name
    char numbuf[16];
    sprintf(numbuf, "%d", corners );
    msq_std::string name = tagName + numbuf;

      // Try to get existing handle
    cornerHandles[corners] = mesh->tag_get( name, err );  MSQ_ERRZERO(err);
      // If got handle, make sure type is correct
    if (cornerHandles[corners])
    {
      Mesh::TagType type;
      size_t size;
      mesh->tag_properties( cornerHandles[corners], name, type, size, err );
      MSQ_ERRZERO(err); 

      if (type != tagType || size != corners*tagLen)
      {
        MSQ_SETERR(err)("Invalid tag type", MsqError::INVALID_ARG );
        return 0;
      }
    }
      // If didn't get handle, try to create it
    else 
    {
      cornerHandles[corners] = mesh->tag_create( name, tagType, tagLen * corners, 0, err );
      MSQ_ERRZERO(err);
    }
  }
  
  return cornerHandles[corners];
}

void CornerTagHandles::save_load_tags( bool load, PatchData* pd, void* data, size_t bytes, MsqError& err )
{
  unsigned i;
  Mesh* mesh = get_current_mesh( pd );
    
    // We want to group the reads into as many elements at
    // once as we can.  We cannot do all the elements at 
    // once because different tags must be used for different 
    // element types.  So group by blocks of the same type of
    // element.  Typically patches contain only one element type
    // so most of the time this loop will be run only once.
  i = 0;
  MsqMeshEntity* elem_array = pd->get_element_array( err ); MSQ_ERRRTN(err);
  while (i < pd->num_elements())
  {
    unsigned start = i++;
    unsigned num_corners = elem_array[start].vertex_count();
    while (i < pd->num_elements() && elem_array[i].vertex_count() == num_corners)
      ++i;
      // elements in [start,i) are the same type
    
      // Find the tag handle
    TagHandle tag = mesh ? get_handle( mesh, num_corners, err ) : 0; MSQ_ERRRTN(err);
   
      // Read/write the data
    const unsigned count = i - start;
    const unsigned offset = pd->get_element_corner_offset(start);
    char* byte_ptr = reinterpret_cast<char*>(data) + offset * bytes;
    const Mesh::ElementHandle* handles = pd->get_element_handles_array() + start;

      // Need to support local cached values w/out a mesh to store
      // the tags on.  If there is no mesh, skip the save/load step.
    if (mesh)
    {
      if (load)
        mesh->tag_get_element_data( tag, count, handles, byte_ptr, err );
      else
        mesh->tag_set_element_data( tag, count, handles, byte_ptr, err );
      MSQ_ERRRTN(err);
    }
  }
}

} // namespace Mesquite
