/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MeshDomain1D.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MeshDomain1D.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

    
void PointDomain::snap_to( Mesh::VertexHandle,
                           Vector3D &coordinate) const
  { coordinate = geom(); }
  
void PointDomain::vertex_normal_at( Mesh::VertexHandle ,
                                    Vector3D &coordinate) const
  { coordinate.set(0,0,0); }

void PointDomain::element_normal_at( Mesh::ElementHandle ,
                                     Vector3D &coordinate) const
  { coordinate.set(0,0,0); }
              
void PointDomain::vertex_normal_at( const Mesh::VertexHandle* ,
                                    Vector3D coordinates[],
                                    unsigned count,
                                    MsqError& err ) const
{ 
  std::fill( coordinates, coordinates+count, Vector3D(0,0,0) );
  MSQ_SETERR(err)( "Cannot get normal for PointDomain", MsqError::INTERNAL_ERROR );
}
                
void PointDomain::closest_point( Mesh::VertexHandle ,
                                 const Vector3D& ,
                                 Vector3D& closest,
                                 Vector3D& normal,
                                 MsqError& err ) const
{
  closest = geom();
  normal.set(0,0,0);
  MSQ_SETERR(err)( "Cannot get normal for PointDomain", MsqError::INTERNAL_ERROR );
}

void PointDomain::domain_DoF( const Mesh::VertexHandle* ,
                              unsigned short* dof_array,
                              size_t num_handles,
                              MsqError&  ) const
  { std::fill( dof_array, dof_array+num_handles, 0 ); }
    
void LineDomain::snap_to( Mesh::VertexHandle ,
                          Vector3D &coordinate) const
  { coordinate = geom().point(geom().closest( coordinate )); }

void LineDomain::vertex_normal_at( Mesh::VertexHandle ,
                                   Vector3D &coordinate) const
    // no normal, return tangent instead.
  { coordinate = geom().direction(); }

void LineDomain::element_normal_at( Mesh::ElementHandle ,
                                    Vector3D &coordinate) const
    // no normal, return tangent instead.
  { coordinate = geom().direction(); }
              
void LineDomain::vertex_normal_at( const Mesh::VertexHandle* ,
                                   Vector3D coordinates[],
                                   unsigned count,
                                   MsqError& err ) const
{
  std::fill( coordinates, coordinates+count, geom().direction() );
  MSQ_SETERR(err)( "Cannot get normal for LineDomain", MsqError::INTERNAL_ERROR );
}
                
void LineDomain::closest_point( Mesh::VertexHandle ,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const
{
  closest = geom().point(geom().closest( position ));
  normal = geom().direction();
  MSQ_SETERR(err)( "Cannot get normal for LineDomain", MsqError::INTERNAL_ERROR );
}
            
                    
void LineDomain::domain_DoF( const Mesh::VertexHandle* ,
                             unsigned short* dof_array,
                             size_t num_handles,
                             MsqError& err ) const
  { std::fill( dof_array, dof_array+num_handles, 1 ); }

void CircleDomain::snap_to( Mesh::VertexHandle ,
                            Vector3D &coordinate) const
  { coordinate = geom().closest( coordinate ); }

void CircleDomain::vertex_normal_at( Mesh::VertexHandle ,
                                     Vector3D &coordinate) const
{
    // no normal, return tangent instead.
  Vector3D junk, copy(coordinate);
  if (!geom().closest( copy, junk, coordinate )) // at center?
    coordinate = geom().radial_vector();
}
  

void CircleDomain::element_normal_at( Mesh::ElementHandle h,
                                      Vector3D &coordinate) const
  { CircleDomain::vertex_normal_at( h, coordinate ); }
  
void CircleDomain::vertex_normal_at( const Mesh::VertexHandle* handles,
                                     Vector3D coordinates[],
                                     unsigned count,
                                     MsqError& err ) const
{
  for (unsigned i = 0; i < count; ++i)
    vertex_normal_at( handles[i], coordinates[i] );
  MSQ_SETERR(err)( "Cannot get normal for CircleDomain", MsqError::INTERNAL_ERROR );
}
  
void CircleDomain::closest_point( Mesh::VertexHandle ,
                                  const Vector3D& position,
                                  Vector3D& closest,
                                  Vector3D& normal,
                                  MsqError& err ) const
{
    // no normal, get tangent instead
  if (!geom().closest(position, closest, normal)) { // at center?
    normal = geom().radial_vector();
    closest = geom().center() + normal;
  }
  MSQ_SETERR(err)( "Cannot get normal for CircleDomain", MsqError::INTERNAL_ERROR );
}
  
void CircleDomain::domain_DoF( const Mesh::VertexHandle* ,
                               unsigned short* dof_array,
                               size_t num_handles,
                               MsqError& err ) const
  { std::fill( dof_array, dof_array+num_handles, 1 ); }



} // namespace Mesquite
