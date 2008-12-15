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

    kraftche@cae.wisc.edu    

  ***************************************************************** */
/** \file QuadLagrangeShape.cpp
 *  \author Jason Kraftcheck
 */
 
#include "QuadLagrangeShape.hpp"
#include "MsqError.hpp"

namespace Mesquite {

EntityTopology QuadLagrangeShape::element_topology() const
  { return QUADRILATERAL; }

void QuadLagrangeShape::coefficients( unsigned loc_dim,
                                      unsigned loc_num,
                                      unsigned nodebits,
                                      double* coeff_out,
                                      size_t& num_coeff,
                                      MsqError& err ) const
{
  num_coeff = 9;
  switch (loc_dim) {
    case 0:
      coeff_out[0] = coeff_out[1] = coeff_out[2] =
      coeff_out[3] = coeff_out[4] = coeff_out[5] = 
      coeff_out[6] = coeff_out[7] = coeff_out[8] = 0.0;
      coeff_out[loc_num] = 1;
      break;
    case 1:
      coeff_out[0] = coeff_out[1] = coeff_out[2] =
      coeff_out[3] = coeff_out[4] = coeff_out[5] = 
      coeff_out[6] = coeff_out[7] = coeff_out[8] = 0.0;
      if (nodebits & (1 << loc_num)) {  
          // if mid-edge node is present
        coeff_out[loc_num+4] = 1;
      }
      else {
          // If mid-edge node is not present, mapping function value
          // for linear edge is even weight of adjacent vertices.
        coeff_out[ loc_num     ] = 0.5;
        coeff_out[(loc_num+1)%4] = 0.5;
      }
      break;
    case 2:
      if (nodebits & 1<<4) { // if center node is present
        coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] =
        coeff_out[4] = coeff_out[5] = coeff_out[6] = coeff_out[7] = 0.0;
        coeff_out[8] = 1;
      } 
      else {
          // for linear element, (no mid-edge nodes), all corners contribute 1/4.
        coeff_out[0] = 0.25;
        coeff_out[1] = 0.25;
        coeff_out[2] = 0.25;
        coeff_out[3] = 0.25;
        coeff_out[8] = 0.00;
          // add in contribution for any mid-edge nodes present
        for (int i = 0; i < 4; ++i) { // for each edge
          if (nodebits & (1<<i))
          {
            coeff_out[ i+4   ]  = 0.5;
            coeff_out[ i     ] -= 0.25;
            coeff_out[(i+1)%4] -= 0.25;
          }
          else {
            coeff_out[ i+4   ]  = 0.0;
          }
        }
      }
      break;
    default:
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                  "Request for dimension %d mapping function value"
                  "for a quadrilateral element", loc_dim);
  }
}
     

static void derivatives_at_corner( unsigned corner, 
                                   unsigned nodebits,
                                   size_t* vertices,
                                   MsqVector<2>* derivs,
                                   size_t& num_vtx )
{
  static const unsigned xi_adj_corners[]  = { 1, 0, 3, 2 };
  static const unsigned xi_adj_edges[]    = { 0, 0, 2, 2 };
  static const unsigned eta_adj_corners[] = { 3, 2, 1, 0 };
  static const unsigned eta_adj_edges[]   = { 3, 1, 1, 3 };
  
  static const double corner_xi[]  = { -0.5,  0.5,  0.5, -0.5 }; // xi values by corner
  static const double corner_eta[] = { -0.5, -0.5,  0.5,  0.5 }; // eta values by corner
  static const double other_xi[]   = {  0.5, -0.5, -0.5,  0.5 }; // xi values for adjacent corner in xi direction
  static const double other_eta[]  = {  0.5,  0.5, -0.5, -0.5 }; // eta values for adjcent corner in eta direction
  static const double mid_xi[]     = {  1.0, -1.0, -1.0,  1.0 }; // xi values for mid-node in xi direction
  static const double mid_eta[]    = {  1.0,  1.0, -1.0, -1.0 }; // eta values for mid-node in eta direction
  
  num_vtx = 3;
  vertices[0] = corner;
  vertices[1] = xi_adj_corners[corner];
  vertices[2] = eta_adj_corners[corner];

  derivs[0][0] = corner_xi [corner];
  derivs[0][1] = corner_eta[corner];
  derivs[1][0] = other_xi  [corner];
  derivs[1][1] = 0.0;
  derivs[2][0] = 0.0;
  derivs[2][1] = other_eta [corner];

  if (nodebits & (1<<xi_adj_edges[corner])) {
    vertices[num_vtx] = 4 + xi_adj_edges[corner];
    derivs[num_vtx][0] = 2.0*mid_xi[corner];
    derivs[num_vtx][1] = 0.0;
    derivs[0][0] -= mid_xi[corner];
    derivs[1][0] -= mid_xi[corner];
    ++num_vtx;
  }

  if (nodebits & (1<<eta_adj_edges[corner])) {
    vertices[num_vtx] = 4 + eta_adj_edges[corner];
    derivs[num_vtx][0] = 0.0;
    derivs[num_vtx][1] = 2.0*mid_eta[corner];
    derivs[0][1] -= mid_eta[corner];
    derivs[2][1] -= mid_eta[corner];
    ++num_vtx;
  }
}

static void derivatives_at_mid_edge( unsigned edge, 
                                     unsigned nodebits,
                                     size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
  static const double values[][9] = { {-0.25, -0.25, 0.25,  0.25, -0.5,  1.0,  0.5,  1.0,  2.0},
                                      {-0.25,  0.25, 0.25, -0.25, -1.0,  0.5, -1.0, -0.5, -2.0},
                                      {-0.25, -0.25, 0.25,  0.25, -0.5, -1.0,  0.5, -1.0, -2.0},
                                      {-0.25,  0.25, 0.25, -0.25,  1.0,  0.5,  1.0, -0.5,  2.0} };
  static const double edge_values[][2] = { {-0.5,  0.5},
                                           {-0.5,  0.5},
                                           { 0.5, -0.5},
                                           { 0.5, -0.5} }; 
  const unsigned prev_corner = edge;           // index of start vertex of edge
  const unsigned next_corner = (edge+1)%4;     // index of end vertex of edge
  const unsigned is_eta_edge = edge % 2;       // edge is xi = +/- 0
  const unsigned is_xi_edge  = 1 - is_eta_edge;// edge is eta = +/- 0
  //const unsigned mid_edge_node = edge + 4;     // mid-edge node index
  const unsigned prev_opposite = (prev_corner+3)%4; // index of corner adjacent to prev_corner
  const unsigned next_opposite = (next_corner+1)%4; // index of corner adjacent to next_corner;
 
    // First do derivatives along edge (e.g. wrt xi if edge is eta = +/-1)
  num_vtx = 2;
  vertices[0] = prev_corner;
  vertices[1] = next_corner;
  derivs[0][is_eta_edge] = edge_values[edge][0];
  derivs[0][is_xi_edge]  = 0.0;
  derivs[1][is_eta_edge] = edge_values[edge][1];
  derivs[1][is_xi_edge]  = 0.0;
    // That's it for the edge-direction derivatives.  No other vertices contribute.
    
    // Next handle the linear element case.  Handle this as a special case first,
    // so the generalized solution doesn't impact performance for linear elements
    // too much.
  if (!nodebits) {
    num_vtx = 4;
    vertices[2] = prev_opposite;
    vertices[3] = next_opposite;
    derivs[0][is_xi_edge] = values[edge][prev_corner];
    derivs[1][is_xi_edge] = values[edge][next_corner];
    derivs[2][is_xi_edge] = values[edge][prev_opposite];
    derivs[2][is_eta_edge] = 0.0;
    derivs[3][is_xi_edge] = values[edge][next_opposite];
    derivs[3][is_eta_edge] = 0.0;
    return;
  }
  
    // Initial (linear) contribution for each corner
  double v[4] = { values[edge][0], 
                  values[edge][1], 
                  values[edge][2], 
                  values[edge][3] };

    // If mid-face node is present
  double v8 = 0.0;
  if (nodebits & 16u) {
    v8 = values[edge][8];
    vertices[num_vtx] = 8;
    derivs[num_vtx][is_eta_edge] = 0.0;
    derivs[num_vtx][is_xi_edge] = v8;
    v[0] -= 0.25 * v8;
    v[1] -= 0.25 * v8;
    v[2] -= 0.25 * v8;
    v[3] -= 0.25 * v8;
    ++num_vtx;
  }

    // If mid-edge nodes are present
  for (unsigned i = 0; i < 4; ++i) {
    if (nodebits & (1<<i)) {
      const double value = values[edge][i+4] - 0.5 * v8;
      if (fabs(value) > 0.125) {
        v[ i     ] -= 0.5 * value;
        v[(i+1)%4] -= 0.5 * value;
        vertices[num_vtx] = i+4;
        derivs[num_vtx][is_eta_edge] = 0.0;
        derivs[num_vtx][is_xi_edge] = value;
        ++num_vtx;
      }
    }
  }

    // update values for adjacent corners
  derivs[0][is_xi_edge] = v[prev_corner];
  derivs[1][is_xi_edge] = v[next_corner];
    // do other two corners
  if (fabs(v[prev_opposite]) > 0.125) {
    vertices[num_vtx] = prev_opposite;
    derivs[num_vtx][is_eta_edge] = 0.0;
    derivs[num_vtx][is_xi_edge] = v[prev_opposite];
    ++num_vtx;
  }
  if (fabs(v[next_opposite]) > 0.125) {
    vertices[num_vtx] = next_opposite;
    derivs[num_vtx][is_eta_edge] = 0.0;
    derivs[num_vtx][is_xi_edge] = v[next_opposite];
    ++num_vtx;
  }
}


static void derivatives_at_mid_elem( unsigned nodebits,
                                     size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
    // fast linear case
    // This is provided as an optimization for linear elements.
    // If this block of code were removed, the general-case code
    // below should produce the same result.
  if (!nodebits) {
    num_vtx = 4;
    vertices[0] = 0; derivs[0][0] = -0.25; derivs[0][1] = -0.25;
    vertices[1] = 1; derivs[1][0] =  0.25; derivs[1][1] = -0.25;
    vertices[2] = 2; derivs[2][0] =  0.25; derivs[2][1] =  0.25;
    vertices[3] = 3; derivs[3][0] = -0.25; derivs[3][1] =  0.25;
    return;
  }
  
  const unsigned n4bit = 1<<0;
  const unsigned n5bit = 1<<1;
  const unsigned n6bit = 1<<2;
  const unsigned n7bit = 1<<3;
  
  num_vtx = 0;
  
    // N_0
  if ((nodebits&(n4bit|n7bit)) != (n4bit|n7bit)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 0;
    derivs[num_vtx][0] = (nodebits&n7bit) ? 0.0 : -0.25;
    derivs[num_vtx][1] = (nodebits&n4bit) ? 0.0 : -0.25;
    ++num_vtx;
  }
  
    // N_1
  if ((nodebits&(n4bit|n5bit)) != (n4bit|n5bit)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 1;
    derivs[num_vtx][0] = (nodebits&n5bit) ? 0.0 :  0.25;
    derivs[num_vtx][1] = (nodebits&n4bit) ? 0.0 : -0.25;
    ++num_vtx;
  }
  
    // N_2
  if ((nodebits&(n5bit|n6bit)) != (n5bit|n6bit)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 2;
    derivs[num_vtx][0] = (nodebits&n5bit) ? 0.0 :  0.25;
    derivs[num_vtx][1] = (nodebits&n6bit) ? 0.0 :  0.25;
    ++num_vtx;
  }
  
    // N_3
  if ((nodebits&(n6bit|n7bit)) != (n6bit|n7bit)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 3;
    derivs[num_vtx][0] = (nodebits&n7bit) ? 0.0 : -0.25;
    derivs[num_vtx][1] = (nodebits&n6bit) ? 0.0 :  0.25;
    ++num_vtx;
  }
  
    // N_4
  if (nodebits&n4bit) {
    vertices[num_vtx] = 4;
    derivs[num_vtx][0] =  0.0;
    derivs[num_vtx][1] = -0.5;
    ++num_vtx;
  }
  
    // N_5
  if (nodebits&n5bit) {
    vertices[num_vtx] = 5;
    derivs[num_vtx][0] =  0.5;
    derivs[num_vtx][1] =  0.0;
    ++num_vtx;
  }
  
    // N_6
  if (nodebits&n6bit) {
    vertices[num_vtx] = 6;
    derivs[num_vtx][0] =  0.0;
    derivs[num_vtx][1] =  0.5;
    ++num_vtx;
  }
  
    // N_7
  if (nodebits&n7bit) {
    vertices[num_vtx] = 7;
    derivs[num_vtx][0] = -0.5;
    derivs[num_vtx][1] =  0.0;
    ++num_vtx;
  }
  
    // N_8 (mid-quad node) never contributes to Jacobian at element center!!!
}

void QuadLagrangeShape::derivatives( unsigned loc_dim,
                                     unsigned loc_num,
                                     unsigned nodebits,
                                     size_t* vertex_indices_out,
                                     MsqVector<2>* d_coeff_d_xi_out,
                                     size_t& num_vtx,
                                     MsqError& err ) const
{
  switch (loc_dim) {
    case 0:
      derivatives_at_corner( loc_num, nodebits, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 1:
      derivatives_at_mid_edge( loc_num, nodebits, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 2:
      derivatives_at_mid_elem( nodebits, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}


} // namespace Mesquite
