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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MappingFunction.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MappingFunction.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

namespace Mesquite {

void MappingFunction2D::jacobian( const PatchData& pd,
                                  size_t element_number,
                                  unsigned nodebits,
                                  unsigned loc_dim,
                                  unsigned loc_num,
                                  size_t* vertex_patch_indices_out,
                                  MsqVector<2>* d_coeff_d_xi_out,
                                  size_t& num_vtx_out,
                                  MsqMatrix<3,2>& jacobian_out,
                                  MsqError& err ) const
{
  const MsqMeshEntity& elem = pd.element_by_index( element_number );
  const size_t* conn = elem.get_vertex_index_array();
  
  derivatives( loc_dim, loc_num, nodebits, vertex_patch_indices_out,
               d_coeff_d_xi_out, num_vtx_out, err ); MSQ_ERRRTN(err);
 
  PatchData::reduced_connectivity_map( elem.get_element_type(),
                                       elem.node_count(),
                                       num_vtx_out,
                                       vertex_patch_indices_out,
                                       vertex_patch_indices_out,
                                       err );  MSQ_ERRRTN(err);
 
  jacobian_out.zero();
  size_t w = 0;
  for (size_t r = 0; r < num_vtx_out; ++r) {
    size_t i = vertex_patch_indices_out[r] = conn[vertex_patch_indices_out[r]];
    MsqMatrix<3,1> coords( pd.vertex_by_index( i ).to_array() );
    jacobian_out += coords * transpose(d_coeff_d_xi_out[r]);
    if (i < pd.num_free_vertices()) {
      vertex_patch_indices_out[w] = i;
      d_coeff_d_xi_out[w] = d_coeff_d_xi_out[r];
      ++w;
    }
  }
  num_vtx_out = w;
}

void MappingFunction3D::jacobian( const PatchData& pd,
                                  size_t element_number,
                                  unsigned nodebits,
                                  unsigned loc_dim,
                                  unsigned loc_num,
                                  size_t* vertex_patch_indices_out,
                                  MsqVector<3>* d_coeff_d_xi_out,
                                  size_t& num_vtx_out,
                                  MsqMatrix<3,3>& jacobian_out,
                                  MsqError& err ) const
{
  const MsqMeshEntity& elem = pd.element_by_index( element_number );
  const size_t* conn = elem.get_vertex_index_array();
  
  derivatives( loc_dim, loc_num, nodebits, vertex_patch_indices_out,
               d_coeff_d_xi_out, num_vtx_out, err ); MSQ_ERRRTN(err);
 
  PatchData::reduced_connectivity_map( elem.get_element_type(),
                                       elem.node_count(),
                                       num_vtx_out,
                                       vertex_patch_indices_out,
                                       vertex_patch_indices_out,
                                       err );  MSQ_ERRRTN(err);
 
  jacobian_out.zero();
  size_t w = 0;
  for (size_t r = 0; r < num_vtx_out; ++r) {
    size_t i = vertex_patch_indices_out[r] = conn[vertex_patch_indices_out[r]];
    MsqMatrix<3,1> coords( pd.vertex_by_index( i ).to_array() );
    jacobian_out += coords * transpose(d_coeff_d_xi_out[r]);
    if (i < pd.num_free_vertices()) {
      vertex_patch_indices_out[w] = i;
      d_coeff_d_xi_out[w] = d_coeff_d_xi_out[r];
      ++w;
    }
  }
  num_vtx_out = w;
}


} // namespace Mesquite
