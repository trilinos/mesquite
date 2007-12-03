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


/** \file data.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "data.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"

#include <iostream>
#include <fstream>

using namespace Mesquite;

const char* temp_file = "syncrononous_input.vtk";
void create_input_mesh( double mid_x , MeshImpl& mesh, MsqError& err )
{
  std::ofstream vtkfile( temp_file );
  vtkfile << "# vtk DataFile Version 3.0" << std::endl
          << "Mesquite Syncronous Boundary test" << std::endl
          << "ASCII" << std::endl
          << "DATASET UNSTRUCTURED_GRID" << std::endl
          << "POINTS 9 float" << std::endl
          << min_x << ' ' << max_y << ' ' << z << std::endl
          << mid_x << ' ' << max_y << ' ' << z << std::endl
          << max_x << ' ' << max_y << ' ' << z << std::endl
          << min_x << ' ' << mid_y << ' ' << z << std::endl
          << mid_x << ' ' << mid_y << ' ' << z << std::endl
          << max_x << ' ' << mid_y << ' ' << z << std::endl
          << min_x << ' ' << min_y << ' ' << z << std::endl
          << mid_x << ' ' << min_y << ' ' << z << std::endl
          << max_x << ' ' << min_y << ' ' << z << std::endl
          << "CELLS 4 20" << std::endl
          << "4 1 0 3 4" << std::endl
          << "4 2 1 4 5" << std::endl
          << "4 4 3 6 7" << std::endl
          << "4 5 4 7 8" << std::endl
          << "CELL_TYPES 4" << std::endl
          << "9 9 9 9" << std::endl
          << "POINT_DATA 9" << std::endl
          << "SCALARS fixed int" << std::endl
          << "LOOKUP_TABLE default" << std::endl
          << "1 0 1" << std::endl
          << "1 0 1" << std::endl
          << "1 0 1" << std::endl
          ;
          
  mesh.read_vtk( temp_file, err );
  remove( temp_file );
}

void MyDomain::setup( Mesh* mesh, MsqError& err )
{
  mHandles.clear();
  mHandles.resize(9);
  
  std::vector<Mesh::EntityHandle> vertices;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRRTN(err);
  if (vertices.size() != 9) {
    MSQ_SETERR(err)(MsqError::INVALID_MESH, "Expected 9 vertices, got %d", (int)vertices.size() );
    return;
  }
  
  MsqVertex coords[9];
  mesh->vertices_get_coordinates( &vertices[0], coords, 9, err ); MSQ_ERRRTN(err);
  
    // figure out which vertices are which
  const double epsilon = 1e-3;
  for (int i = 0; i < 9; ++i) {
    if (fabs(coords[i][2] - z) > epsilon) {
      MSQ_SETERR(err)(MsqError::INVALID_MESH, "Invalid Z coord" );
      return;
    }
    
    int x, y;
    if (fabs(coords[i][1] - min_y) < epsilon)
      y = 2;
    else if (fabs(coords[i][1] - mid_y) < epsilon)
      y = 1;
    else if (fabs(coords[i][1] - max_y) < epsilon)
      y = 0;
    else {      
      MSQ_SETERR(err)(MsqError::INVALID_MESH, "Invalid Y coord" );
      return;
    }
    
    if (fabs(coords[i][0] - min_x) < epsilon)
      x = 0;
    else if (fabs(coords[i][0] - max_x) < epsilon)
      x = 2;
    else if (min_x < coords[i][0] && max_x > coords[i][0])
      x = 1;
    else {      
      MSQ_SETERR(err)(MsqError::INVALID_MESH, "Invalid X coord: %f", coords[i][0] );
      return;
    }
    
    mHandles[3*y+x] = vertices[i];
  }
}

void MyDomain::snap_to(Mesh::EntityHandle entity_handle,
                       Vector3D &coordinate) const
{
    // everything gets moved into the plane
  coordinate[2] = z;
    // for boundary vertices, check that they are on appropriate curve/vertex
  switch (index(entity_handle)) {
    case 0: coordinate[0] = min_x; coordinate[1] = max_y; break;
    case 1:                        coordinate[1] = max_y; break;
    case 2: coordinate[0] = max_x; coordinate[1] = max_y; break;

    case 3: coordinate[0] = min_x;                        break;
    case 5: coordinate[0] = max_x;                        break;

    case 6: coordinate[0] = min_x; coordinate[1] = min_y; break;
    case 7:                        coordinate[1] = min_y; break;
    case 8: coordinate[0] = max_x; coordinate[1] = min_y; break;
  }
}
  
void MyDomain::normal_at(Mesh::EntityHandle handle, Vector3D &norm) const
{
  norm.set(0,0,1);
}

void MyDomain::normal_at( const Mesh::EntityHandle* handles,
                          Vector3D normals[],
                          unsigned count,
                          MsqError&  ) const
{
  for (unsigned i = 0; i < count; ++i)
    normal_at( handles[i], normals[i] );
}
    
void MyDomain::closest_point( Mesh::EntityHandle handle,
                              const Vector3D& position,
                              Vector3D& closest,
                              Vector3D& normal,
                              MsqError&  ) const
{
  normal = position;
  normal_at( handle, normal );
  closest = position;
  snap_to( handle, closest );
}
    
void MyDomain::domain_DoF( const Mesh::EntityHandle* handle_array,
                           unsigned short* dof_array,
                           size_t num_handles,
                           MsqError&  ) const
{
  for (unsigned i = 0; i < num_handles; ++i) {
    switch (index(handle_array[i])) {
        // corners
      case 0: 
      case 2: 
      case 6:
      case 8:
        dof_array[i] = 0;
        break;
        
        //edges:
      case 1: 
      case 3: 
      case 5:
      case 7:
        dof_array[i] = 1;
        break;
      
        // every thing else (middle vertex, quads)
      default:
        dof_array[i] = 2;
        break;
    }
  }
}

