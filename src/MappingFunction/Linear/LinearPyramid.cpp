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


/** \file LinearPyramid.cpp
 *  \brief LinearPyramid implementation
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LinearPyramid.hpp"
#include "MsqError.hpp"

namespace Mesquite {

static const char* nonlinear_error 
 = "Attempt to use LinearTriangle mapping function for a nonlinear element\n";

static inline void set_equal_derivatives( double value, 
                                          size_t* indices,
                                          MsqVector<3>* derivs,
                                          size_t& num_vtx )
{
  num_vtx = 5;
  indices[0] = 0;
  indices[1] = 1;
  indices[2] = 2;
  indices[3] = 3;
  indices[4] = 4;
    
  derivs[0][0] = -value;
  derivs[0][1] = -value;
  derivs[0][2] = -0.125;
  
  derivs[1][0] =  value;
  derivs[1][1] = -value;
  derivs[1][2] = -0.125;
  
  derivs[2][0] =  value;
  derivs[2][1] =  value;
  derivs[2][2] = -0.125;
  
  derivs[3][0] = -value;
  derivs[3][1] =  value;
  derivs[3][2] = -0.125;
  
  derivs[4][0] =  0.0;
  derivs[4][1] =  0.0;
  derivs[4][2] =  0.5;
}

static inline void set_edge_derivatives( unsigned base_corner,
                                         double value,
                                         size_t* indices,
                                         MsqVector<3>* derivs,
                                         size_t& num_vtx )
{
  const int direction = base_corner % 2;
  const int edge_beg =  base_corner;
  const int edge_end = (base_corner+1)%4;
  const int adj_end  = (base_corner+2)%4;
  const int adj_beg  = (base_corner+3)%4;
  const int dir_sign = 2*(edge_beg/2) - 1;
  const int oth_sign = 2*((edge_beg+1)/2%2) - 1;

  num_vtx = 5;
  indices[0] = edge_beg;
  indices[1] = edge_end;
  indices[2] = adj_end;
  indices[3] = adj_beg;
  indices[4] = 4;

  derivs[0][  direction] =  2 * dir_sign * value;
  derivs[0][1-direction] =      oth_sign * value;
  derivs[0][2] = -0.25;

  derivs[1][  direction] = -2 * dir_sign * value;
  derivs[1][1-direction] =      oth_sign * value;
  derivs[1][2] = -0.25;

  derivs[2][  direction] =  0.0;
  derivs[2][1-direction] = -oth_sign * value;
  derivs[2][2]           =  0.0;

  derivs[3][  direction] =  0.0;
  derivs[3][1-direction] = -oth_sign * value;
  derivs[3][2]           =  0.0;

  derivs[4][0] = 0.0;
  derivs[4][1] = 0.0;
  derivs[4][2] = 0.5;
}

static inline void set_corner_derivatives( unsigned corner,
                                           double value,
                                           size_t* indices,
                                           MsqVector<3>* derivs,
                                           size_t& num_vtx )
{
  const unsigned adj_in_xi = (5 - corner) % 4;
  const unsigned adj_in_eta = 3 - corner;

  const int dxi_sign  = 2*((corner+1)/2%2)-1;
  const int deta_sign = 2*(corner/2) - 1;
  const double dxi_value = dxi_sign * value;
  const double deta_value = deta_sign * value;

  num_vtx = 4;
  indices[0] = corner;
  indices[1] = adj_in_xi;
  indices[2] = adj_in_eta;
  indices[3] = 4;

  derivs[0][0] =  dxi_value;
  derivs[0][1] =  deta_value;
  derivs[0][2] = -0.5;

  derivs[1][0] = -dxi_value;
  derivs[1][1] =  0.0;
  derivs[1][2] =  0.0;

  derivs[2][0] =  0.0;
  derivs[2][1] = -deta_value;
  derivs[2][2] =  0.0;

  derivs[3][0] =  0.0;
  derivs[3][1] =  0.0;
  derivs[3][2] =  0.5;
}

EntityTopology LinearPyramid::element_topology() const
  { return PYRAMID; }
  
int LinearPyramid::num_nodes() const
  { return 5; }

static void coefficients_at_corner( unsigned corner,
                                    double* coeff_out,
                                    size_t* indices_out,
                                    size_t& num_coeff )
{
  num_coeff = 1;
  indices_out[0] = corner;
  coeff_out[0] = 1.0;
}


static void coefficients_at_mid_edge( unsigned edge,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  num_coeff = 2;
  coeff_out[0] = 0.5;
  coeff_out[1] = 0.5;
  
  if (edge < 4) {
    indices_out[0] = edge;
    indices_out[1] = (edge+1)%4;
  }
  else {
    indices_out[0] = edge-4;
    indices_out[1] = 4;
  }
}

static void coefficients_at_mid_face( unsigned face,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  if (face == 4) {
    num_coeff = 4;
    coeff_out[0] = 0.25;
    coeff_out[1] = 0.25;
    coeff_out[2] = 0.25;
    coeff_out[3] = 0.25;
    indices_out[0] = 0;
    indices_out[1] = 1;
    indices_out[2] = 2;
    indices_out[3] = 3;
  }
  else {
    num_coeff = 3;
    indices_out[0] = face;
    indices_out[1] = (face+1)%4;
    indices_out[2] = 4;
    coeff_out[0] = MSQ_ONE_THIRD;
    coeff_out[1] = MSQ_ONE_THIRD;
    coeff_out[2] = MSQ_ONE_THIRD;
  }
}

static void coefficients_at_mid_elem( double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  num_coeff = 5;
  coeff_out[0] = 0.1875;
  coeff_out[1] = 0.1875;
  coeff_out[2] = 0.1875;
  coeff_out[3] = 0.1875;
  coeff_out[4] = 0.25;
  indices_out[0] = 0;
  indices_out[1] = 1;
  indices_out[2] = 2;
  indices_out[3] = 3;
  indices_out[4] = 4;
}

void LinearPyramid::coefficients( unsigned loc_dim,
                                  unsigned loc_num,
                                  NodeSet nodeset,
                                  double* coeff_out,
                                  size_t* indices_out,
                                  size_t& num_coeff,
                                  MsqError& err ) const
{
  if (nodeset.have_any_mid_node()) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  switch (loc_dim) {
    case 0:
      coefficients_at_corner( loc_num, coeff_out, indices_out, num_coeff );
      break;
    case 1:
      coefficients_at_mid_edge( loc_num, coeff_out, indices_out, num_coeff );
      break;
    case 2:
      coefficients_at_mid_face( loc_num, coeff_out, indices_out, num_coeff );
      break;
    case 3:
      coefficients_at_mid_elem( coeff_out, indices_out, num_coeff );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

void LinearPyramid::derivatives( unsigned loc_dim,
                                 unsigned loc_num,
                                 NodeSet nodeset,
                                 size_t* vertex_indices_out,
                                 MsqVector<3>* d_coeff_d_xi_out,
                                 size_t& num_vtx,
                                 MsqError& err ) const
{
  if (nodeset.have_any_mid_node()) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  switch (loc_dim) {
    case 0:
      if (loc_num == 4) {
        set_equal_derivatives( 0.0, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      else {
        set_corner_derivatives( loc_num, 0.5, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      break;
    case 1:
      if (loc_num < 4) {
        set_edge_derivatives( loc_num, 0.25, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      else {
        set_corner_derivatives( loc_num-4, 0.25, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }    
      break;
    case 2:
      if (loc_num == 4) {
        set_equal_derivatives( 0.25, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      else {
        set_edge_derivatives( loc_num, 1./6, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      break;
    case 3:
      set_equal_derivatives( 0.1875, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

} // namespace Mesquite
