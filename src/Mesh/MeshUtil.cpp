/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MeshUtil.cpp
 *  \brief Implement MeshUtil class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MeshUtil.hpp"
#include "MeshInterface.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"

#include <vector>
#include <algorithm>
#include <limits>

namespace MESQUITE_NS {


void MeshUtil::edge_length_distribution( double& min_out,
                                         double& avg_out,
                                         double& rms_out,
                                         double& max_out,
                                         double& std_dev_out,
                                         MsqError& err )
{
  size_t count = 0;
  avg_out = rms_out = std_dev_out = 0.0;
  min_out = std::numeric_limits<double>::max();
  max_out = std::numeric_limits<double>::min();

    // some lists 
  std::vector<Mesh::VertexHandle> verts, conn, adj;
  std::vector<Mesh::VertexHandle>::iterator i;
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::VertexHandle>::iterator j;
  std::vector<size_t> junk(2);
  std::vector<MsqVertex> coords;

    // get all vertices in the mesh
  myMesh->get_all_vertices( verts, err ); MSQ_ERRRTN(err);
  
    // for each vertex in the mesh
  for (i = verts.begin(); i != verts.end(); ++i) {
      // initially clear what will eventually be the list of 
      // all opposite ends of edges at vertex i
    adj.clear();
  
      // get all adjacent elements
    elems.clear();
    myMesh->vertices_get_attached_elements( &*i, 1, elems, junk, err ); MSQ_ERRRTN(err);
    
      // for each adjacent element
    for (j = elems.begin(); j != elems.end(); ++j) {
      EntityTopology type;
      myMesh->elements_get_topologies( &*j, &type, 1, err ); MSQ_ERRRTN(err);
      conn.clear();
      myMesh->elements_get_attached_vertices( &*j, 1, conn, junk, err ); MSQ_ERRRTN(err);
  
        // get position of input vertex in connectivity list
      size_t idx = std::find( conn.begin(), conn.end(), *i ) - conn.begin();
      if (idx > TopologyInfo::corners(type))
        continue; // skip higher-order nodes
  
        // get the set of other corner vertices connected to vertex i
        // by an edge
      unsigned n;
      const unsigned* indices = TopologyInfo::adjacent_vertices( type, idx, n );
      
        // keep only those vertex handles for vertices that are
        // a) connected to this vertex by an edge and b) for which the
        // handle value for this vertex is less than the other, adjacent
        // vertex (this avoids counting edges more than once.)
      for (unsigned k = 0; k < n; ++k)
        if (*i < conn[indices[k]])
          adj.push_back( conn[indices[k]] );
    }
    
      // remove duplicates
    std::sort( adj.begin(), adj.end() );
    adj.erase( std::unique( adj.begin(), adj.end() ), adj.end() );
    adj.push_back( *i );
    
      // now calculate edge lengths
    coords.resize( adj.size() );
    myMesh->vertices_get_coordinates( &adj[0], &coords[0], adj.size(), err ); MSQ_ERRRTN(err);
    for (size_t i = 0; i < coords.size()-1; ++i) {
      Vector3D diff = coords[i] - coords.back();
      double len_sqr = diff % diff;
      double len = sqrt(len_sqr);
      avg_out += len;
      rms_out += len_sqr;
      if (len < min_out)
        min_out = len;
      if (len > max_out)
        max_out = len;
      ++count;
    }
  }
  
  if (!count) { // no mesh
    MSQ_SETERR(err)("Mesh contains no elements", MsqError::INVALID_MESH);
    return;
  }

    // calculate final stats
  avg_out /= count;
  rms_out /= count;
  std_dev_out = sqrt( rms_out - avg_out*avg_out );
  rms_out = sqrt(rms_out);
}



} // namespace MESQUITE_NS
