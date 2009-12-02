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
  \file   MeshTransform.cpp
  \brief  

  The MeshTransform Class is the base class for all the smoothing algorythms 

  \author Michael Brewer
  \date   2004-11-06
*/


#include "MeshTransform.hpp"
#include "MeshInterface.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {
/*! Constructor
\param in_mat Matrix component of transform.  Specifies rotation,
              scaling, and reflection.
\param in_vec Vector component of transform.  Specifies translation.
*/
  MeshTransform::MeshTransform(Matrix3D &in_mat, Vector3D &in_vec)
  {
    mMat = in_mat;
    mVec = in_vec;
  }
  
/*! 
  Actually apply the affine transformation
  */
  double MeshTransform::loop_over_mesh( Mesh* mesh,
                                        MeshDomain* ,
                                        const Settings*,
                                        MsqError &err )
  {
    std::vector<Mesh::VertexHandle> handle_list;
    mesh->get_all_vertices( handle_list, err );
    if (MSQ_CHKERR(err))
      return 1.0;
    
    MsqVertex vertex;
    std::vector<Mesh::VertexHandle>::const_iterator iter;
    for (iter = handle_list.begin(); iter != handle_list.end(); ++iter)
    {
      mesh->vertices_get_coordinates( &*iter, &vertex, 1, err );
      if (MSQ_CHKERR(err))
        return 1.0;
      
      vertex = mMat * vertex + mVec;
      
      mesh->vertex_set_coordinates( *iter, vertex, err );
      if (MSQ_CHKERR(err))
        return 1.0;
    }

    return 0.0;
  }
  
  
  void MeshTransform::add_translation( const Vector3D& offset )
    { mVec += offset; }
  
  void MeshTransform::add_rotation( const Vector3D& a, double radians )
    {
      const double c = cos(radians);
      const double s = sin(radians);
      const Matrix3D m1(    c,   -a[2]*s, a[1]*s,
                          a[2]*s,   c,   -a[0]*s,
                         -a[1]*s, a[0]*s,   c    );
      mMat = m1 * mMat;
      mVec = m1 * mVec;
    }
  
  void MeshTransform::add_scale( double factor )
    { add_scale( Vector3D(factor) ); } 
  
  void MeshTransform::add_scale( const Vector3D& f )
    {
      for (int i = 0; i < 3; ++i) {
        mVec[i] *= f[i];
        mMat[i][0] *= f[i];
        mMat[i][1] *= f[i];
        mMat[i][2] *= f[i];
      }
    }      
  
} // namespace Mesquite
