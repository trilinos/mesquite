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

#ifndef MSQ_CYLINDER_DOMAIN_CPP
#define MSQ_CYLINDER_DOMAIN_CPP

#include "Mesquite.hpp"
#include "CylinderDomain.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <limits.h>
#else
#  include <limits>
#endif

namespace Mesquite {

void CylinderDomain::evaluate( Mesh::EntityHandle,
                               const Vector3D& point,
                               Vector3D& closest,
                               Vector3D& normal ) const
{
  const double EPSILON = msq_std::numeric_limits<double>::epsilon();
  double t = mAxis % (point - mCenter);
  const Vector3D axis_point = mCenter + t * mAxis;
  
  normal = point - axis_point;
  const double len = normal.length();
  if (len < EPSILON)
  {
    Vector3D v( 1, 0, 0 );
    if ((v * mAxis).length() < EPSILON)
      v.set( 0, 1, 0 );
    normal = v * mAxis;
    normal.normalize();
  }
  else
  {
    normal /= len;
  }
  
  closest = axis_point + mRadius * normal;
}

void CylinderDomain::snap_to( Mesh::EntityHandle h, Vector3D& v ) const
{
  Vector3D p(v), n;
  evaluate( h, p, v, n );
}

void CylinderDomain::normal_at( Mesh::EntityHandle h, Vector3D& v ) const
{
  Vector3D p(v), l;
  evaluate( h, p, l, v );
}

void CylinderDomain::normal_at( const Mesh::EntityHandle* h,
                                Vector3D* coords,
                                unsigned count, 
                                MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i)
    normal_at( h[i], coords[i] );
}
 
void CylinderDomain::closest_point( Mesh::EntityHandle handle,
                                    const Vector3D& position,
                                    Vector3D& closest,
                                    Vector3D& normal,
                                    MsqError&  ) const
{
  evaluate( handle, position, closest, normal );
}

void CylinderDomain::domain_DoF( const Mesh::EntityHandle* ,
                                 unsigned short* dof_array,
                                 size_t count,
                                 MsqError&  ) const
{
  msq_std::fill( dof_array, dof_array + count, (unsigned short)2 );
}
  


} // namespace Mesquite

#endif
