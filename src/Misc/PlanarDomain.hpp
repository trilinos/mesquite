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
/*!
  \file   PlanarDomain.hpp
  \brief  


  \author Thomas Leurent
  \date   2002-01-17
*/


#ifndef MSQ_PLANAR_DOMAIN_HPP
#define MSQ_PLANAR_DOMAIN_HPP

#include "MeshInterface.hpp"

namespace Mesquite
{
  class Mesh;
  
  /*! \class PlanarDomain
       This is a template for a planar domain.
       It will provide the normal information necessary for surface mesh optimization.
    */
  class PlanarDomain : public Mesquite::MeshDomain
  {
  public:
    PlanarDomain(const Vector3D& normal,
                 const Vector3D& point,
                 Mesquite::Mesh* mesh)
        : mNormal(normal / normal.length()),
          mPoint(point),
          mMesh(mesh)
      {}
    
    virtual ~PlanarDomain() { }

    void set_plane(const Vector3D& normal, const Vector3D& point)
      {
        mNormal = normal / normal.length();
        mPoint = point;
      }

    void set_mesh(Mesquite::Mesh* mesh)
      { mMesh = mesh; }
    
    virtual void snap_to(Mesh::EntityHandle entity_handle,
                         Vector3D &coordinate);
    
    virtual void normal_at(Mesh::EntityHandle entity_handle,
                           Vector3D &coordinate);

  private:
    Vector3D mNormal;
    Vector3D mPoint;
    Mesquite::Mesh* mMesh;
  };
}

#endif
