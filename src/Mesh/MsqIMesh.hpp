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


/** \file MsqIMesh.hpp
 *  \brief Adaptor for ITAPS iMesh interface
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_I_MESH_HPP
#define MSQ_I_MESH_HPP

#include "MeshInterface.hpp"
#include "iMesh.h"

namespace Mesquite {

/** The name of the tag (integer) that Mesquite will use
 *  to store internal data
 */
const char* const VERTEX_BYTE_TAG_NAME  = "MesquiteVertexByte";

/** The name of the tag (integer) Mesquite expects to be non-zero
 *  for vertices which are not to be moved by Mesquite
 */
const char* const VERTEX_FIXED_TAG_NAME = "MesquiteVertexFixed";

/** The name of the tag (integer) Mesquite expects to be non-zero
 *  for vertices that are higher-order nodes slaved to their logical
 *  position.
 */
const char* const VERTEX_SLAVED_TAG_NAME = "MesquiteVertexSlaved";

/**\class MsqIMesh
 *\brief Mesquite iMesh Adapter
 *
 * Adpater for interfacing Mesquite with an application that provides
 * the ITAPS iMesh interface for interacting with mesh data.
 */
class MsqIMesh : public Mesquite::Mesh
{
public:

  /**\brief factory method
   *
   * Create an MsqIMesh instance.
   */
  static MsqIMesh* create( iMesh_Instance imesh, 
                           iBase_EntitySetHandle meshset,
                           iBase_EntityType element_dimension,
                           MsqError& err,
                           const char* fixed_tag_name = VERTEX_FIXED_TAG_NAME,
                           const char* slaved_tag_name= VERTEX_SLAVED_TAG_NAME );
  
  /**\brief factory method
   *
   * Create an MsqIMesh instance with no initial mesh.
   * Call set_active_set to complete initialization.
   */
  static MsqIMesh* create( iMesh_Instance imesh, 
                           MsqError& err,
                           const char* fixed_tag_name = VERTEX_FIXED_TAG_NAME,
                           const char* slaved_tag_name= VERTEX_SLAVED_TAG_NAME);
  
  virtual void set_active_set( iBase_EntitySetHandle meshset, 
                               iBase_EntityType element_dimension,
                               MsqError& err ) = 0;
  
  virtual iMesh_Instance get_imesh_instance() const = 0;
  virtual iBase_EntitySetHandle get_entity_set() const = 0;
  
  virtual ~MsqIMesh() = 0;
};  


} // namespace Mesquite

#endif
