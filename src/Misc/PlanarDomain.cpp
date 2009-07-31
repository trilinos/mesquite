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
#include "PlanarDomain.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <algorithm.h>
#else
# include <algorithm>
#endif

void Mesquite::PlanarDomain::set_plane( const Mesquite::Vector3D& normal, 
                                        const Mesquite::Vector3D& point)
{
  mNormal = normal;
  mNormal.normalize();
  mCoeff = -(mNormal % point);
}

void Mesquite::PlanarDomain::snap_to(Mesquite::Mesh::VertexHandle
                                       entity_handle,
                                     Vector3D &coordinate) const
{
  coordinate -= mNormal * ( mNormal % coordinate + mCoeff );
}


void Mesquite::PlanarDomain::vertex_normal_at(
  Mesquite::Mesh::VertexHandle /*entity_handle*/,
  Mesquite::Vector3D &coordinate) const
{
  coordinate = mNormal;
}

void Mesquite::PlanarDomain::element_normal_at(
  Mesquite::Mesh::ElementHandle /*entity_handle*/,
  Mesquite::Vector3D &coordinate) const
{
  coordinate = mNormal;
}


void Mesquite::PlanarDomain::vertex_normal_at( 
                                        const Mesquite::Mesh::VertexHandle* ,
                                        Vector3D coords[],
                                        unsigned count,
                                        Mesquite::MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i)
    coords[i] = mNormal;
}

void Mesquite::PlanarDomain::closest_point( Mesquite::Mesh::VertexHandle ,
                                            const Mesquite::Vector3D& position,
                                            Mesquite::Vector3D& closest,
                                            Mesquite::Vector3D& normal,
                                            Mesquite::MsqError& ) const
{
  normal = mNormal;
  closest = position - mNormal * (mNormal % position + mCoeff);
}

void Mesquite::PlanarDomain::domain_DoF( const Mesh::VertexHandle* ,
                                         unsigned short* dof_array,
                                         size_t num_vertices,
                                         MsqError&  ) const
{
  msq_std::fill( dof_array, dof_array + num_vertices, 2 );
}

Mesquite::PlanarDomain Mesquite::PlanarDomain::fit_vertices( Mesquite::Mesh* mesh, 
                                                             double epsilon, 
                                                             Mesquite::MsqError& err )
{
  // get all vertex coordiantes
  msq_std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err ); 
  if (MSQ_CHKERR(err))
    return PlanarDomain(XY);
  if (verts.size() < 3) {
    MSQ_SETERR(err)("No mesh", MsqError::INVALID_MESH);
    return PlanarDomain(XY);
  }
  msq_std::vector<MsqVertex> coords(verts.size());
  mesh->vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err ); 
  if (MSQ_CHKERR(err))
    return PlanarDomain(XY);
  
    // Assume all vertices are co-planer.  Find a triple of
    // vertices that will result in a small rounding error.
  Vector3D box_min(HUGE_VAL,HUGE_VAL,HUGE_VAL),
           box_max(-HUGE_VAL,-HUGE_VAL,-HUGE_VAL);
  for (size_t i = 0; i < coords.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (box_min[j] > coords[i][j])
        box_min[j] = coords[i][j];
      if (box_max[j] < coords[i][j])
        box_max[j] = coords[i][j];
    }
  }
  Vector3D corners[8] = { Vector3D(box_min[0],box_min[1],box_min[2]),
                          Vector3D(box_max[0],box_min[1],box_min[2]),
                          Vector3D(box_min[0],box_max[1],box_min[2]),
                          Vector3D(box_max[0],box_max[1],box_min[2]),
                          Vector3D(box_min[0],box_min[1],box_max[2]),
                          Vector3D(box_max[0],box_min[1],box_max[2]),
                          Vector3D(box_min[0],box_max[1],box_max[2]),
                          Vector3D(box_max[0],box_max[1],box_max[2]) };
                          
  assert(sizeof(char) == sizeof(bool));
  std::vector<unsigned char> fixed( verts.size() );
  mesh->vertices_get_fixed_flag( &verts[0], (bool*)&fixed[0], verts.size(), err );
  size_t first = std::find(fixed.begin(), fixed.end(), 1) - fixed.begin();
    // Our goal here is to consider only the boundary (fixed) vertices
    // when calculating the plane.  If for some reason the user wants
    // to snap a not-quite-planar mesh to a plane during optimization, 
    // if possible we want to start with the plane containing the fixed
    // vertices, as those won't get snapped.  If no vertices are fixed,
    // then treat them all as fixed for the purpose calculating the plane
    // (consider them all.)
  if (first == fixed.size()) {
    first = 0;
    std::fill( fixed.begin(), fixed.end(), 1 );
  }
  
  Vector3D pts[8] = { coords[first], coords[first], coords[first], coords[first],
                      coords[first], coords[first], coords[first], coords[first] };
  for (size_t i = first+1; i < coords.size(); ++i) {
    for (int j = 0; j < 8; ++j) {
      if (fixed[i] && 
          ((pts[j] - corners[j]).length_squared() > 
           (coords[i] - corners[j]).length_squared()))
        pts[j] = coords[i];
    }
  }
  
  Vector3D normal(0,0,0);
  for (int i = 0; i < 6; ++i) 
    for (int j = i+1; j < 7; ++j)
      for (int k = j+1; k < 8; ++k) {
        Vector3D n = (pts[j]-pts[i]) * (pts[k]-pts[i]);
        if (n.length_squared() > normal.length_squared())
          normal = n;
      }
  
  if (normal.length_squared() < epsilon*epsilon) {
    MSQ_SETERR(err)("All vertices colinear", MsqError::INVALID_MESH);
    return PlanarDomain(XY);
  }
  
    // check that all vertices line in plane
  Vector3D point = coords[0];
  normal /= normal.length();
  const double d = -(normal % point);
  for (size_t i = 0; i < coords.size(); ++i) {
    double dist = fabs( normal % coords[i] + d );
    if (dist > epsilon) {
      MSQ_SETERR(err)("Vertices are not coplanar", MsqError::INVALID_MESH );
      return PlanarDomain(XY);
    }
  }
  
    // free memory used to hold vertex data
  msq_std::vector<MsqVertex>* tmp_vvect = new msq_std::vector<MsqVertex>;
  tmp_vvect->swap( coords );
  delete tmp_vvect;
  msq_std::vector<Mesh::VertexHandle>* tmp_hvect = new msq_std::vector<Mesh::VertexHandle>;
  tmp_hvect->swap( verts );
  delete tmp_hvect;
  
    // now count inverted elements
  coords.resize(3);
  size_t inverted_count = 0;
  msq_std::vector<Mesh::ElementHandle> elems;
  msq_std::vector<size_t> junk;
  mesh->get_all_elements( elems, err ); 
  if (MSQ_CHKERR(err))
    return PlanarDomain(XY);
  for (size_t i = 0; i < elems.size(); ++i) {
    
    verts.clear();
    mesh->elements_get_attached_vertices( &elems[0], 1, verts, junk, err );
    if (MSQ_CHKERR(err))
      return PlanarDomain(XY);

    EntityTopology type;
    mesh->elements_get_topologies( &elems[i], &type, 1, err );
    if (TopologyInfo::dimension(type) != 2 || verts.size() < 3) {
      MSQ_SETERR(err)("Mesh contains non-surface elements", MsqError::INVALID_MESH);
      return PlanarDomain(XY);
    }
    
    mesh->vertices_get_coordinates( &verts[0], &coords[0], 3, err ); 
    if (MSQ_CHKERR(err))
      return PlanarDomain(XY);
    Vector3D n = (coords[1] - coords[0]) * (coords[2] - coords[0]);
    if (n % normal < 0.0)
      ++inverted_count;
  }
  
    // if most elements are inverted, flip normal
  if (2*inverted_count > elems.size())
    normal = -normal;
  
  return PlanarDomain( normal, point );
}
    

    
