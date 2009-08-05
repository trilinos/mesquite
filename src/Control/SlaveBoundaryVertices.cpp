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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file SlaveBoundaryVertices.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "SlaveBoundaryVertices.hpp"
#include "Settings.hpp"
#include "MsqError.hpp"
#include "MeshInterface.hpp"
#include "MsqVertex.hpp"

#include <vector>
#include <algorithm>

namespace MESQUITE_NS {

SlaveBoundaryVertices::SlaveBoundaryVertices( unsigned depth, unsigned dim )
  : elemDepth(depth), domainDoF(dim)
  {}

msq_std::string SlaveBoundaryVertices::get_name() const
  { return "SlaveBoundaryVertices"; }

struct BoolArr { // vector<bool> cannot be treated as an array of bool
  BoolArr(size_t size) : mArray(new bool[size]) {}
  ~BoolArr() { delete [] mArray; }
  bool* mArray;
  bool& operator[](size_t i) { return mArray[i]; }
};

double SlaveBoundaryVertices::loop_over_mesh( Mesh* mesh, 
                                              MeshDomain* domain, 
                                              const Settings* settings,
                                              MsqError& err )
{
  if (settings->get_slaved_ho_node_mode() != Settings::SLAVE_CALCULATED) {
    MSQ_SETERR(err)("Request to calculate higher-order node slaved status "
                    "when Settings::get_get_slaved_ho_node_mode() "
                    "!= SLAVE_CALCUALTED", MsqError::INVALID_STATE);
    return 0.0;
  }
  
    // If user said that we should treat fixed vertices as the
    // boundary, but specified that fixed vertices are defined
    // by the dimension of their geometric domain, then just
    // do distance from the geometric domain.
  int dim = this->domainDoF;
  if (dim >= 4 && settings->get_fixed_vertex_mode() != Settings::FIXED_FLAG)
    dim = settings->get_fixed_vertex_mode();
  
    // Create a map to contain vertex depth.  Intiliaze all to 
    // elemDepth+1.
  msq_std::vector<Mesh::VertexHandle> vertices;
  mesh->get_all_vertices( vertices, err );  MSQ_ERRZERO(err);
  if (vertices.empty())
    return 0.0;
  msq_std::sort( vertices.begin(), vertices.end() );
  msq_std::vector<unsigned short> depth( vertices.size(), elemDepth+1 );
  BoolArr fixed( vertices.size() );
  mesh->vertices_get_fixed_flag( &vertices[0], fixed.mArray, vertices.size(), err );
  MSQ_ERRZERO(err);
  
    // Initialize map with boundary vertices.
  if (dim >= 4) {
    for(size_t i = 0; i < vertices.size(); ++i)
      if (fixed[i])
        depth[i] = 0;
  }
  else {
    if (!domain) {
      MSQ_SETERR(err)("Request to calculate higher-order node slaved status "
                      "by distance from bounding domain without a domain.",
                      MsqError::INVALID_STATE);
      return 0.0;
    }
    
    msq_std::vector<unsigned short> dof( vertices.size() );
    domain->domain_DoF( &vertices[0], &dof[0], vertices.size(), err ); MSQ_ERRZERO(err);
    for (size_t i = 0; i < vertices.size(); ++i)
      if (dof[i] <= dim)
        depth[i] = 0;
  }
  
    // Now iterate over elements repeatedly until we've found all of the 
    // elements near the boundary.  This could be done much more efficiently
    // using vertex-to-element adjacencies, but it is to common for 
    // applications not to implement that unless doing relaxation smoothing.
    // This is O(elemDepth * elements.size() * ln(vertices.size()));
  msq_std::vector<Mesh::ElementHandle> elements;
  msq_std::vector<Mesh::ElementHandle>::const_iterator j, k;
  msq_std::vector<Mesh::VertexHandle> conn;
  msq_std::vector<size_t> junk(2);
  mesh->get_all_elements( elements, err );  MSQ_ERRZERO(err);
  if (elements.empty())
    return 0.0;
  bool some_changed;
  do {
    some_changed = false;
    for (j = elements.begin(); j != elements.end(); ++j) {
      conn.clear(); junk.clear();
      mesh->elements_get_attached_vertices( &*j, 1, conn, junk, err ); MSQ_ERRZERO(err);
      unsigned short elem_depth = elemDepth+1;
      for (k = conn.begin(); k != conn.end(); ++k) {
        size_t i = msq_std::lower_bound( vertices.begin(), vertices.end(), *k ) - vertices.begin();
        if (i == vertices.size()) {
          MSQ_SETERR(err)("Invalid vertex handle in element connectivity list.", 
                          MsqError::INVALID_MESH);
          return 0.0;
        }
        if (depth[i] < elem_depth)
          elem_depth = depth[i];
      }
      if (elem_depth == elemDepth+1)
        continue;
      
      ++elem_depth;
      for (k = conn.begin(); k != conn.end(); ++k) {
        size_t i = msq_std::lower_bound( vertices.begin(), vertices.end(), *k ) - vertices.begin();
        if (depth[i] > elem_depth) {
          depth[i] = elem_depth;
          some_changed = true;
        }
      }
    } // for(elements)
  } while (some_changed);
  
    // Now remove any corner vertices from the slaved set
  msq_std::vector<Mesh::VertexHandle>::iterator p;
  msq_std::vector<EntityTopology> types(elements.size());
  mesh->elements_get_topologies( &elements[0], &types[0], elements.size(), err ); MSQ_ERRZERO(err);
  for (j = elements.begin(); j != elements.end(); ++j) {
    const unsigned corners = TopologyInfo::corners(types[j-elements.begin()]);
    conn.clear(); junk.clear();
    mesh->elements_get_attached_vertices( &*j, 1, conn, junk, err ); MSQ_ERRZERO(err);
    for (unsigned i = 0; i < corners; ++i) {
      p = msq_std::lower_bound( vertices.begin(), vertices.end(), conn[i] );
      depth[p-vertices.begin()] = 0;
    }
  }
  
    // Now mark all vertices *not* within specified depth as slave vertices.
  msq_std::vector<unsigned char> bytes( vertices.size() );
  mesh->vertices_get_byte( &vertices[0], &bytes[0], vertices.size(), err ); MSQ_ERRZERO(err);
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (depth[i] <= elemDepth || fixed[i])
      bytes[i] &= ~MsqVertex::MSQ_DEPENDENT;
    else
      bytes[i] |= MsqVertex::MSQ_DEPENDENT;
  }
  
  mesh->vertices_set_byte( &vertices[0], &bytes[0], vertices.size(), err ); MSQ_ERRZERO(err);
  return 0.0;
}


} // namespace Mesquite
