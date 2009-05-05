/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
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

/*!
  \file   MsqIRel.cpp
  \brief  


  \author Jason Kraftcheck
  \date   2007-08-14
*/

#include "MsqIRel.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "MsqIBase.hpp"

namespace MESQUITE_NS
{


/***************** DomainTSTT class methods *********************/

MsqIRel::MsqIRel( iGeom_Instance geom,
                  iRel_Instance relate_iface,
                  iRel_RelationHandle relate_instance ) 
  : MsqCommonIGeom( geom ), 
    relateIface( relate_iface ),
    relateInstance( relate_instance )
{
}

MsqIRel::~MsqIRel() {}


void MsqIRel::snap_to( Mesh::VertexHandle handle,
                           Vector3D& coordinate ) const
{
  int ierr;
  iBase_EntityHandle geom;
  
  ierr = geom_from_mesh( handle, geom );
  if (iBase_SUCCESS != ierr) {
    process_itaps_error( ierr );
    return;
  }
  
  ierr = move_to( geom, coordinate );
  if (iBase_SUCCESS != ierr) {
    process_itaps_error( ierr );
    return;
  }
}

void MsqIRel::vertex_normal_at( Mesh::VertexHandle handle,
                                    Vector3D& coordinate ) const
{
  int ierr;
  iBase_EntityHandle geom;
  
  ierr = geom_from_mesh( handle, geom );
  if (iBase_SUCCESS != ierr) {
    process_itaps_error( ierr );
    return;
  }
  
  ierr = normal( geom, coordinate );
  if (iBase_SUCCESS != ierr) {
    process_itaps_error( ierr );
    return;
  }
}

void MsqIRel::element_normal_at( Mesh::ElementHandle handle,
                                     Vector3D& coordinate ) const
{
  MsqIRel::vertex_normal_at( handle, coordinate );
}

void MsqIRel::vertex_normal_at( const Mesh::VertexHandle* handle,
                                    Vector3D coordinates[],
                                    unsigned count,
                                    MsqError& err ) const
{
  int ierr;
  
  geomHandles.resize( count );
  ierr = geom_from_mesh( handle, &geomHandles[0], count );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
  
  ierr = normal( &geomHandles[0], coordinates, count );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
}

void MsqIRel::domain_DoF( const Mesh::VertexHandle* handle_array,
                              unsigned short* dof_array,
                              size_t count,
                              MsqError& err ) const
{
  int ierr;
  
  geomHandles.resize( count );
  ierr = geom_from_mesh( handle_array, &geomHandles[0], count );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
  
  ierr = get_dimension( &geomHandles[0], dof_array, count );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
}
    


void MsqIRel::closest_point( Mesh::VertexHandle handle,
                                 const Vector3D& position,
                                 Vector3D& closest,
                                 Vector3D& normal,
                                 MsqError& err ) const
{
  int ierr;
  iBase_EntityHandle geom;
  
  ierr = geom_from_mesh( handle, geom );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
  
  ierr = closest_and_normal( geom, position, closest, normal );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
}


int MsqIRel::geom_from_mesh( Mesh::EntityHandle mesh_ent_handle,
                                 iBase_EntityHandle& geom_handle ) const
{
    // get geometric entity
  int ierr;
  iRel_getEntEntAssociation( relateIface,
                             relateInstance,
                             (iBase_EntityHandle)mesh_ent_handle,
                             false,
                             &geom_handle,
                             &ierr );
  if (iBase_SUCCESS != ierr)
    return ierr;
  
    // get dimension of geometric entities
  int type, one = 1, one_too = 1, *type_ptr = &type;
  iGeom_getArrType( geomIFace, &geom_handle, 1, &type_ptr, &one, &one_too, &ierr );
  if (iBase_SUCCESS != ierr)
    return ierr;
  
    // not interested in volumes (only surfaces, curves, and points)
  if (type == iBase_REGION)
    geom_handle = 0;
  
  return iBase_SUCCESS;
}


int MsqIRel::geom_from_mesh( const Mesh::EntityHandle* handles,
                                 iBase_EntityHandle* geom_handles,
                                 size_t count ) const
{
  int ierr;
  for (size_t i = 0; i < count; ++i) {
    ierr = geom_from_mesh( handles[i], geom_handles[i] );
    if (iBase_SUCCESS != ierr)
      return ierr;
  }
  
  return iBase_SUCCESS;
}                    


} // namespace Mesquite

