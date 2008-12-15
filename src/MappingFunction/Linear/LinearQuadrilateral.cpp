/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */
/** \file LinearQuadrilateral.cpp
 *  \author Jason Kraftcheck
 */
 
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "LinearQuadrilateral.hpp"

namespace Mesquite {

static const char* nonlinear_error 
 = "Attempt to use LinearQuadrilateral mapping function for a nonlinear element\n";

EntityTopology LinearQuadrilateral::element_topology() const
  { return QUADRILATERAL; }

void LinearQuadrilateral::coefficients( unsigned loc_dim,
                                        unsigned loc_num,
                                        unsigned nodebits,
                                        double* coeff_out,
                                        size_t& num_coeff,
                                        MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 4;
  switch (loc_dim) {
    case 0:
      coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] = 0.0;
      coeff_out[loc_num] = 1.0;
      break;
    case 1:
      coeff_out[ loc_num     ] = coeff_out[(loc_num+1)%4] = 0.5;
      coeff_out[(loc_num+2)%4] = coeff_out[(loc_num+3)%4] = 0.0;
      break;
    case 2:
      coeff_out[0] = 0.25;
      coeff_out[1] = 0.25;
      coeff_out[2] = 0.25;
      coeff_out[3] = 0.25;
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

const unsigned  xi = 0;
const unsigned eta = 1;
const int sign[2][4] = {{ -1,  1,  1, -1 },  // xi
                        { -1, -1,  1,  1 }}; // eta

static void derivatives_at_corner( unsigned corner, 
                                   size_t* vertex_indices,
                                   MsqVector<2>* derivs,
                                   size_t& num_vtx )
{
  const unsigned adj_in_xi = (5 - corner) % 4;
  const unsigned adj_in_eta = 3 - corner;
  
  num_vtx = 3;
  vertex_indices[0] = corner;
  vertex_indices[1] = adj_in_xi;
  vertex_indices[2] = adj_in_eta;
  
  derivs[0][0] = 0.5 * sign[ xi][corner];
  derivs[0][1] = 0.5 * sign[eta][corner];
  derivs[1][0] = 0.5 * sign[ xi][adj_in_xi ];
  derivs[1][1] = 0.0;
  derivs[2][0] = 0.0;
  derivs[2][1] = 0.5 * sign[eta][adj_in_eta];
}

static void derivatives_at_mid_edge( unsigned edge, 
                                     size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
  const unsigned start_vtx =  edge;
  const unsigned   end_vtx = (edge+1) % 4;
  const unsigned othr1_vtx = (edge+2) % 4;
  const unsigned othr2_vtx = (edge+3) % 4;
  const unsigned direction = edge % 2;
  const unsigned orthogonal = 1 - direction;
  
  num_vtx = 4;
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 2;
  vertices[3] = 3;
  
  derivs[start_vtx][direction] = 0.5 * sign[direction][start_vtx];
  derivs[  end_vtx][direction] = 0.5 * sign[direction][  end_vtx];
  derivs[othr1_vtx][direction] = 0.0;
  derivs[othr2_vtx][direction] = 0.0;
 
  derivs[0][orthogonal] = 0.25 * sign[orthogonal][0];
  derivs[1][orthogonal] = 0.25 * sign[orthogonal][1];
  derivs[2][orthogonal] = 0.25 * sign[orthogonal][2];
  derivs[3][orthogonal] = 0.25 * sign[orthogonal][3];
}


static void derivatives_at_mid_elem( size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
  num_vtx = 4;
  vertices[0] = 0; derivs[0][0] = -0.25; derivs[0][1] = -0.25;
  vertices[1] = 1; derivs[1][0] =  0.25; derivs[1][1] = -0.25;
  vertices[2] = 2; derivs[2][0] =  0.25; derivs[2][1] =  0.25;
  vertices[3] = 3; derivs[3][0] = -0.25; derivs[3][1] =  0.25;
}

void LinearQuadrilateral::derivatives( unsigned loc_dim,
                                       unsigned loc_num,
                                       unsigned nodebits,
                                       size_t* vertex_indices_out,
                                       MsqVector<2>* d_coeff_d_xi_out,
                                       size_t& num_vtx,
                                       MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  switch (loc_dim) {
    case 0:
      derivatives_at_corner( loc_num, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 1:
      derivatives_at_mid_edge( loc_num, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 2:
      derivatives_at_mid_elem( vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

} // namespace Mesquite
