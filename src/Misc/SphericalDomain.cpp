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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu         
   
  ***************************************************************** */
#include "Mesquite.hpp"
#include "SphericalDomain.hpp"
#include "Vector3D.hpp"

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <algorithm.h>
#else
# include <algorithm>
#endif

void Mesquite::SphericalDomain::snap_to(Mesh::VertexHandle /*entity_handle*/,
                                        Vector3D &coordinate) const
{
    // Get vector center to coordinate, store in coordinate.
  coordinate -= mCenter;
    // Get distance from center of sphere
  double len = coordinate.length();
    // Scale vector to have length of radius
  coordinate *= mRadius / len;
    // If was at center, return arbitrary position on sphere
    // (all possitions are equally close)
  if (!finite(coordinate.x()))
    coordinate.set( mRadius, 0.0, 0.0 );
    // Get position from vector
  coordinate += mCenter;
}

void Mesquite::SphericalDomain::vertex_normal_at(Mesh::VertexHandle /*entity_handle*/,
                                                 Vector3D &coordinate) const
{
    // normal is vector from center to input position
  coordinate -= mCenter;
    // make it a unit vector
  double length = coordinate.length();
  coordinate /= length;
    // if input position was at center, choose same position
    // on sphere as snap_to.
  if (!finite(coordinate.x()))
    coordinate.set( 1.0, 0.0, 0.0 );
}
void Mesquite::SphericalDomain::element_normal_at(Mesh::ElementHandle h,
                                                 Vector3D &coordinate) const
{
  SphericalDomain::vertex_normal_at( h, coordinate );
}

void Mesquite::SphericalDomain::vertex_normal_at( 
                                const Mesquite::Mesh::VertexHandle* handle,
                                Mesquite::Vector3D coords[],
                                unsigned count,
                                Mesquite::MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i)
    vertex_normal_at( handle[i], coords[i] );
}

void Mesquite::SphericalDomain::closest_point( Mesquite::Mesh::VertexHandle ,
                                               const Mesquite::Vector3D& position,
                                               Mesquite::Vector3D& closest,
                                               Mesquite::Vector3D& normal,
                                               Mesquite::MsqError& ) const
{
  normal = position - mCenter;
  normal.normalize();
  if (!finite(normal.x()))
    normal.set( 1.0, 0.0, 0.0 );
  closest = mCenter + mRadius * normal;
}


void Mesquite::SphericalDomain::domain_DoF( const Mesh::VertexHandle* ,
                                            unsigned short* dof_array,
                                            size_t num_vertices,
                                            MsqError&  ) const
{
  msq_std::fill( dof_array, dof_array + num_vertices, 2 );
}

  
