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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
#include "Mesquite.hpp"
#include "SphericalDomain.hpp"
#include "Vector3D.hpp"

// Currently only works if "entity_handle" refers to a vertex
void Mesquite::SphericalDomain::snap_to(Mesh::EntityHandle /*entity_handle*/,
                                        Vector3D &coordinate)
{
  if (!mMesh)
  {
    return;
  }
  
    //translate sphere to origin
  coordinate -= mCenter;
  
    // Move point to the sphere surface
  double len = coordinate.length();
  
    // If it's not right at the center...
  if (len > MSQ_MIN)
    coordinate *= (mRadius / len);
  
    // Move it back off of the origin
  coordinate += mCenter;
}

void Mesquite::SphericalDomain::normal_at(Mesh::EntityHandle /*entity_handle*/,
                                          Vector3D &coordinate)
{
  if (!mMesh)
  {
    coordinate.set(0, 0, 0);
    return;
  }
  
    //translate sphere to origin
  coordinate -= mCenter;
  
    // See how far it is from the center of the sphere.
  double len=coordinate.length();
  if(len<MSQ_MIN)
  {
      // If it's right at the center...
    coordinate.set(0, 0, 0);
  }
  else
  {
      // Otherwise normalize...
    coordinate /= len;
  }
}
