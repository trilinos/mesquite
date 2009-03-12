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
  \file   MsqIGeom.cpp
  \brief  


  \author Jason Kraftcheck
  \date   2007-08-14
*/

#include "MsqIGeom.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "MsqIBase.hpp"

namespace Mesquite
{

/**\brief Common code for specific implementations of MeshDomain on TSTT interfaces.
 *
 * This class contains the common functionality used by concrete implementations
 * of MeshDomain on the TSTT geometry interface.
 */
class CommonIGeom
{
public:

    /**\param geom The TSTT geometry interface implementation to query */
  CommonIGeom( iGeom_Instance geom );
  
  virtual ~CommonIGeom();

    /** Evaluate the closest point to the input position on the specified
     *  geometric entity and return the result in the passed position 
     *  argument (move the passed position onto the geometry.)
     */
  int move_to( iBase_EntityHandle geom_handle, Vector3D& coord ) const;
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  int normal ( iBase_EntityHandle geom_handle, Vector3D& coord ) const ;
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  int normal ( iBase_EntityHandle geom_handle, Vector3D coords[], unsigned count ) const;
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  int normal ( const iBase_EntityHandle geom_handles[], Vector3D coords[], unsigned count ) const;
    
    /** Given a geometric entity and a position, get point on 
     *  the geometric entity closest to the input position, and
     *  the surface normal at that position.
     */
  int closest_and_normal( iBase_EntityHandle geom_handle,
                           const Vector3D& position,
                           Vector3D& closest, 
                           Vector3D& normal ) const;
                         
  int get_dimension( iBase_EntityHandle const* geom_handle, 
                      unsigned short* dof_out,
                      size_t count ) const ;
                      
  iGeom_Instance geomIFace;

private:  
  mutable std::vector<iBase_EntityHandle> geomHandles;
  mutable std::vector<double> coordArray;
  mutable std::vector<int> typeArray;
};


/* General MeshDomain on iGeom & iRel implementation */

class DomainIGeom : public MsqIGeom, protected CommonIGeom
{
public:

  DomainIGeom( iGeom_Instance geom,
               iRel_Instance irel_iface,
               iRel_RelationHandle irel_instance );

  virtual ~DomainIGeom();

  void snap_to( Mesh::VertexHandle entity_handle,
                Vector3D& coordinat ) const;

  void vertex_normal_at( Mesh::VertexHandle entity_handle,
                         Vector3D& coordinate ) const;

  void element_normal_at( Mesh::ElementHandle entity_handle,
                          Vector3D& coordinate ) const;
  
  void vertex_normal_at( const Mesh::VertexHandle* handles,
                         Vector3D coordinates[],
                         unsigned count,
                         MsqError& err ) const;

  void closest_point( Mesh::VertexHandle handle,
                      const Vector3D& position,
                      Vector3D& closest,
                      Vector3D& normal,
                      MsqError& err ) const;

  void domain_DoF( const Mesh::VertexHandle* handle_array,
                   unsigned short* dof_array,
                   size_t num_vertices,
                   MsqError& err ) const;
                      
protected:

    /** Get geometric entity owning a mesh entity */
  int geom_from_mesh( Mesh::EntityHandle  mesh_handle_in,
                      iBase_EntityHandle& geom_handle_out ) const;
  
  int geom_from_mesh( Mesh::EntityHandle const* mesh_handles_in,
                      iBase_EntityHandle      * geom_handles_out,
                      size_t count ) const;

private:

    /** ITAPS interface implementation for mesh->geometry association */
  iRel_Instance  relateIface;
  iRel_RelationHandle relateInstance;
  
    /** temporary storage of geometry entity handles */
  mutable std::vector<iBase_EntityHandle> geomHandles;
};


/* Specialized one-geometry-entity MeshDomain on iGeom implementation */

class EntIGeom : public MsqIGeom, protected CommonIGeom
{
public:

  EntIGeom( iGeom_Instance geom,
            iBase_EntityHandle geom_ent_handle );

  virtual ~EntIGeom();

  void snap_to( Mesh::VertexHandle entity_handle,
                Vector3D& coordinat ) const;

  void vertex_normal_at( Mesh::VertexHandle entity_handle,
                         Vector3D& coordinate ) const;

  void element_normal_at( Mesh::ElementHandle entity_handle,
                          Vector3D& coordinate ) const;
  
  void vertex_normal_at( const Mesh::VertexHandle* handles,
                         Vector3D coordinates[],
                         unsigned count,
                         MsqError& err ) const;

  void closest_point( Mesh::VertexHandle handle,
                      const Vector3D& position,
                      Vector3D& closest,
                      Vector3D& normal,
                      MsqError& err ) const;

  void domain_DoF( const Mesh::VertexHandle* handle_array,
                   unsigned short* dof_array,
                   size_t num_vertices,
                   MsqError& err ) const;
private:
  
    /** A handle for the geometry entity to evaluate */
  iBase_EntityHandle geomEntHandle;
};




/*****************  GeomTSTT base class methods *********************/

MsqIGeom* MsqIGeom::create( iGeom_Instance geom,
                            iRel_Instance relate_iface,
                            iRel_RelationHandle relate_instance,
                            MsqError& )
{
  return new DomainIGeom( geom, relate_iface, relate_instance );
}

MsqIGeom* MsqIGeom::create( iGeom_Instance geom,
                            iBase_EntityHandle geom_ent_handle,
                            MsqError& )
{
  return new EntIGeom( geom, geom_ent_handle );
}
  

MsqIGeom::~MsqIGeom() {}




/***************** DomainTSTT class methods *********************/

DomainIGeom::DomainIGeom( iGeom_Instance geom,
                          iRel_Instance relate_iface,
                          iRel_RelationHandle relate_instance ) 
  : CommonIGeom( geom ), 
    relateIface( relate_iface ),
    relateInstance( relate_instance )
{
}

DomainIGeom::~DomainIGeom() {}


void DomainIGeom::snap_to( Mesh::VertexHandle handle,
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

void DomainIGeom::vertex_normal_at( Mesh::VertexHandle handle,
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

void DomainIGeom::element_normal_at( Mesh::ElementHandle handle,
                                     Vector3D& coordinate ) const
{
  DomainIGeom::vertex_normal_at( handle, coordinate );
}

void DomainIGeom::vertex_normal_at( const Mesh::VertexHandle* handle,
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

void DomainIGeom::domain_DoF( const Mesh::VertexHandle* handle_array,
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
    


void DomainIGeom::closest_point( Mesh::VertexHandle handle,
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


int DomainIGeom::geom_from_mesh( Mesh::EntityHandle mesh_ent_handle,
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


int DomainIGeom::geom_from_mesh( const Mesh::EntityHandle* handles,
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


/***************** EntIGeom class methods *********************/

EntIGeom::EntIGeom( iGeom_Instance geom, iBase_EntityHandle geom_ent_handle ) 
  : CommonIGeom( geom ), 
    geomEntHandle( geom_ent_handle )
{
}

EntIGeom::~EntIGeom() {}


void EntIGeom::snap_to( Mesh::VertexHandle,
                        Vector3D& coordinate ) const
{
  int ierr = move_to( geomEntHandle, coordinate );
  if (iBase_SUCCESS != ierr)
    process_itaps_error(ierr);
}

void EntIGeom::vertex_normal_at( Mesh::VertexHandle,
                                 Vector3D& coordinate ) const
{
  int ierr = normal( geomEntHandle, coordinate );
  if (iBase_SUCCESS != ierr)
    process_itaps_error(ierr);
}

void EntIGeom::element_normal_at( Mesh::ElementHandle,
                                  Vector3D& coordinate ) const
{
  int ierr = normal( geomEntHandle, coordinate );
  if (iBase_SUCCESS != ierr)
    process_itaps_error(ierr);
}


void EntIGeom::vertex_normal_at( const Mesh::VertexHandle*,
                                 Vector3D coordinates[],
                                 unsigned count,
                                 MsqError& err ) const
{
  int ierr = normal( geomEntHandle, coordinates, count );
  if (iBase_SUCCESS != ierr)
    MSQ_SETERR(err)(process_itaps_error(ierr), MsqError::INTERNAL_ERROR);
}

void EntIGeom::closest_point( Mesh::VertexHandle handle,
                              const Vector3D& position,
                              Vector3D& closest,
                              Vector3D& normal,
                              MsqError& err ) const
{
  int ierr = closest_and_normal( geomEntHandle, position, closest, normal );
  if (iBase_SUCCESS != ierr)
    MSQ_SETERR(err)(process_itaps_error(ierr), MsqError::INTERNAL_ERROR);
}

void EntIGeom::domain_DoF( const Mesh::VertexHandle* ,
                           unsigned short* dof_array,
                           size_t num_vertices,
                           MsqError& err ) const
{
  unsigned short dim;
  int ierr = get_dimension( &geomEntHandle, &dim, 1 );
  if (iBase_SUCCESS != ierr)
    MSQ_SETERR(err)(process_itaps_error(ierr), MsqError::INTERNAL_ERROR);
  msq_std::fill( dof_array, dof_array + num_vertices, dim );
}




/***************** GeomTSTTCommon class methods *********************/

CommonIGeom::CommonIGeom( iGeom_Instance geom )
  : geomIFace( geom )
{
}

CommonIGeom::~CommonIGeom() {}



int CommonIGeom::move_to( iBase_EntityHandle geom, Vector3D& coord ) const
{
  double x, y, z;
  int ierr;
  iGeom_getEntClosestPt( geomIFace, geom, coord[0], coord[1], coord[2], &x, &y, &z, &ierr );
  coord.set( x, y, z );
  return ierr;
}

 
 
int CommonIGeom::normal( iBase_EntityHandle geom, Vector3D& coord ) const
{
  double i, j, k;
  int ierr;
  iGeom_getEntNrmlXYZ( geomIFace, geom, coord[0], coord[1], coord[2], &i, &j, &k, &ierr );
  coord.set( i, j, k );
  return ierr;
}
 
int CommonIGeom::normal( iBase_EntityHandle geom, Vector3D coords[], unsigned count ) const
{
  geomHandles.resize( count, geom );
  return normal( &geomHandles[0], coords, count );
}
 
int CommonIGeom::normal( const iBase_EntityHandle* geom_handles, 
                         Vector3D coords[], 
                         unsigned count ) const
{
    // going to assume this in the following reinterpret_cast, so
    // check to make sure it is true
  assert( sizeof(Vector3D) == 3*sizeof(double) );
  
    // copy input coordinates into array
  coordArray.clear();
  coordArray.resize( 3*count );
  memcpy( &coordArray[0], coords, 3 * count * sizeof(double) );
  
    // define junk variables required for ITAPS "consistancy"
  int junk_1 = count, junk_2 = count;
  double* norm_ptr = reinterpret_cast<double*>(coords);
  
    // get the normals
  int ierr;
  iGeom_getArrNrmlXYZ( geomIFace, 
                       geom_handles,
                       count,
                       iBase_INTERLEAVED,
                       &coordArray[0],
                       count,
                       &norm_ptr,
                       &junk_1,
                       &junk_2,
                       &ierr ); 
  
  return ierr;
}

int CommonIGeom::closest_and_normal( iBase_EntityHandle geom, 
                                     const Vector3D& position,
                                     Vector3D& closest,
                                     Vector3D& normal ) const
{
  int ierr;
  iGeom_getEntNrmlPlXYZ( geomIFace, geom, 
                         position[0], position[1], position[2], 
                         &closest[0], &closest[1], &closest[2],
                          &normal[0],  &normal[1],  &normal[2],
                         &ierr );
  return ierr;
}

                         
int CommonIGeom::get_dimension( const iBase_EntityHandle* geom_handle, 
                                unsigned short* dof_out,
                                size_t count ) const
{
  int ierr;
  typeArray.resize( count );
  
    // define junk variables required for ITAPS "consistancy"
  int junk_1 = count, junk_2 = count;
  int* type_ptr = &typeArray[0];
  
    // get the types
  iGeom_getArrType( geomIFace, geom_handle, count, &type_ptr, &junk_1, &junk_2, &ierr );
  
    // convert from int to unsigned short
  std::copy( typeArray.begin(), typeArray.end(), dof_out );
  return ierr;
}    

} // namespace Mesquite

