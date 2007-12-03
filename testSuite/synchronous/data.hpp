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


/** \file data.hpp
 *  \brief Input data for test
 *  \author Jason Kraftcheck 
 */

#ifndef DATA_HPP
#define DATA_HPP

#include "MeshInterface.hpp"

namespace Mesquite{ class MeshImpl; }

// constants
const double min_x = 0.0, max_x = 2.0;
const double min_y = 0.0, mid_y = 1.0, max_y = 2.0;
const double z = 0.0;

void create_input_mesh( double mid_x, Mesquite::MeshImpl& mesh, Mesquite::MsqError&  );

class MyDomain : public Mesquite::MeshDomain
{
  public:
  
    void setup( Mesquite::Mesh* mesh, Mesquite::MsqError& err );

    void snap_to( Mesquite::Mesh::EntityHandle entity_handle,
                  Mesquite::Vector3D &coordinate) const;
  
    void normal_at( Mesquite::Mesh::EntityHandle entity_handle,
                    Mesquite::Vector3D &coordinate) const;

    void normal_at( const Mesquite::Mesh::EntityHandle* handles,
                    Mesquite::Vector3D coordinates[],
                    unsigned count,
                    Mesquite::MsqError& err ) const;
    
    void closest_point( Mesquite:: Mesh::EntityHandle handle,
                        const Mesquite::Vector3D& position,
                        Mesquite::Vector3D& closest,
                        Mesquite::Vector3D& normal,
                        Mesquite::MsqError& err ) const;
    
    void domain_DoF( const Mesquite::Mesh::EntityHandle* handle_array,
                     unsigned short* dof_array,
                     size_t num_handles,
                     Mesquite::MsqError& err ) const;
 
  private:
  
    std::vector<Mesquite::Mesh::EntityHandle> mHandles;
    inline int index( Mesquite::Mesh::EntityHandle handle ) const;
};

    
int MyDomain::index( Mesquite::Mesh::EntityHandle h ) const
{
  std::vector<Mesquite::Mesh::EntityHandle>::const_iterator i;
  i = std::find( mHandles.begin(), mHandles.end(), h );
  return i == mHandles.end() ? -1 : i - mHandles.begin();
}

#endif
