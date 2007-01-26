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

void QuadLagrangeShape::coefficients_at_corner( unsigned corner,
                                                unsigned ,
                                                msq_std::vector<double>& coeff_out,
                                                MsqError& ) const
{
  coeff_out.clear();
  coeff_out.resize( 9, 0.0 );
  coeff_out[corner] = 1;
}

void QuadLagrangeShape::coefficients_at_mid_edge( unsigned edge,
                                                unsigned nodebits,
                                                msq_std::vector<double>& coeff_out,
                                                MsqError& ) const
{
  coeff_out.clear();
  coeff_out.resize( 9, 0.0 );
  if (nodebits & (1 << edge)) {  
      // if mid-edge node is present
    coeff_out[edge+4] = 1;
  }
  else {
      // If mid-edge node is not present, mapping function value
      // for linear edge is even weight of adjacent vertices.
    coeff_out[ edge     ] = 0.5;
    coeff_out[(edge+1)%4] = 0.5;
  }
}

void QuadLagrangeShape::coefficients_at_mid_face( unsigned ,
                                                unsigned ,
                                                msq_std::vector<double>& ,
                                                MsqError& err) const
{
  MSQ_SETERR(err)("Request for mid-face mapping function value"
                  "for a quadrilateral element", MsqError::UNSUPPORTED_ELEMENT);
}

void QuadLagrangeShape::coefficients_at_mid_elem( unsigned nodebits,
                                          msq_std::vector<double>& coeff_out,
                                          MsqError& ) const
{
  coeff_out.clear();
  coeff_out.resize( 9, 0.0 );
  if (nodebits & 1<<4) { // if center node is present
    coeff_out[8] = 1;
  } 
  else {
      // for linear element, (no mid-edge nodes), all corners contribute 1/4.
    coeff_out[0] = 0.25;
    coeff_out[1] = 0.25;
    coeff_out[2] = 0.25;
    coeff_out[3] = 0.25;
      // add in contribution for any mid-edge nodes present
    for (int i = 0; i < 4; ++i) // for each edge
      if (nodebits & (1<<i))
      {
        coeff_out[ i+4   ]  = 0.5;
        coeff_out[ i     ] -= 0.25;
        coeff_out[(i+1)%4] -= 0.25;
      }
  }
}

void QuadLagrangeShape::derivatives_at_corner( unsigned corner, 
                                               unsigned nodebits,
                                               msq_std::vector<size_t>& vertices,
                                               msq_std::vector<double>& derivs,
                                               MsqError&  ) const
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
  
  vertices.resize(3);
  vertices[0] = corner;
  vertices[1] = xi_adj_corners[corner];
  vertices[2] = eta_adj_corners[corner];

  derivs.resize(6);
  derivs[0] = corner_xi [corner];
  derivs[1] = corner_eta[corner];
  derivs[2] = other_xi  [corner];
  derivs[3] = 0.0;
  derivs[4] = 0.0;
  derivs[5] = other_eta [corner];

  if (nodebits & (1<<xi_adj_edges[corner])) {
    vertices.push_back(4 + xi_adj_edges[corner]);
    derivs.push_back(2.0*mid_xi[corner]);
    derivs.push_back(0.0);
    derivs[0] -= mid_xi[corner];
    derivs[2] -= mid_xi[corner];
  }

  if (nodebits & (1<<eta_adj_edges[corner])) {
    vertices.push_back(4 + eta_adj_edges[corner]);
    derivs.push_back(0.0);
    derivs.push_back(2.0*mid_eta[corner]);
    derivs[1] -= mid_eta[corner];
    derivs[5] -= mid_eta[corner];
  }
}

void QuadLagrangeShape::derivatives_at_mid_edge( unsigned edge, 
                                                 unsigned nodebits,
                                                 msq_std::vector<size_t>& vertices,
                                                 msq_std::vector<double>& derivs,
                                                 MsqError&  ) const
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
  vertices.resize(2);
  vertices[0] = prev_corner;
  vertices[1] = next_corner;
  derivs.clear();
  derivs.resize(4,0.0);
  derivs[0+is_eta_edge] = edge_values[edge][0];
  derivs[2+is_eta_edge] = edge_values[edge][1];
    // That's it for the edge-direction derivatives.  No other vertices contribute.
    
    // Next handle the linear element case.  Handle this as a special case first,
    // so the generalized solution doesn't impact performance for linear elements
    // too much.
  if (!nodebits) {
    vertices.resize(4);
    vertices[2] = prev_opposite;
    vertices[3] = next_opposite;
    derivs.resize(8,0.0);
    derivs[0+is_xi_edge] = values[edge][prev_corner];
    derivs[2+is_xi_edge] = values[edge][next_corner];
    derivs[4+is_xi_edge] = values[edge][prev_opposite];
    derivs[6+is_xi_edge] = values[edge][next_opposite];
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
    vertices.push_back( 8 );
    size_t len = derivs.size();
    derivs.resize( len+2 );
    derivs[len+is_eta_edge] = 0.0;
    derivs[len+is_xi_edge] = v8;
    v[0] -= 0.25 * v8;
    v[1] -= 0.25 * v8;
    v[2] -= 0.25 * v8;
    v[3] -= 0.25 * v8;
  }

    // If mid-edge nodes are present
  for (unsigned i = 0; i < 4; ++i) {
    if (nodebits & (1<<i)) {
      const double value = values[edge][i+4] - 0.5 * v8;
      if (fabs(value) > 0.125) {
        v[ i     ] -= 0.5 * value;
        v[(i+1)%4] -= 0.5 * value;
        vertices.push_back( i+4 );
        size_t len = derivs.size();
        derivs.resize( len+2 );
        derivs[len+is_eta_edge] = 0.0;
        derivs[len+is_xi_edge] = value;
      }
    }
  }

    // update values for adjacent corners
  derivs[0+is_xi_edge] = v[prev_corner];
  derivs[2+is_xi_edge] = v[next_corner];
    // do other two corners
  if (fabs(v[prev_opposite]) > 0.125) {
    vertices.push_back(prev_opposite);
    size_t len = derivs.size();
    derivs.resize(len+2);
    derivs[len+is_eta_edge] = 0.0;
    derivs[len+is_xi_edge] = v[prev_opposite];
  }
  if (fabs(v[next_opposite]) > 0.125) {
    vertices.push_back(next_opposite);
    size_t len = derivs.size();
    derivs.resize(len+2);
    derivs[len+is_eta_edge] = 0.0;
    derivs[len+is_xi_edge] = v[next_opposite];
  }
}

  
void QuadLagrangeShape::derivatives_at_mid_face( unsigned , 
                                                 unsigned ,
                                                 msq_std::vector<size_t>& ,
                                                 msq_std::vector<double>& ,
                                                 MsqError& err ) const

{
  MSQ_SETERR(err)("Request for mid-face mapping function derivative"
                  "for a quadrilateral element", MsqError::UNSUPPORTED_ELEMENT);
}


void QuadLagrangeShape::derivatives_at_mid_elem( unsigned nodebits,
                                                 msq_std::vector<size_t>& vertices,
                                                 msq_std::vector<double>& derivs,
                                                 MsqError& ) const
{
    // fast linear case
    // This is provided as an optimization for linear elements.
    // If this block of code were removed, the general-case code
    // below should produce the same result.
  if (!nodebits) {
    vertices.resize(4);
    derivs.resize(8);
    vertices[0] = 0; derivs[0] = -0.25; derivs[1] = -0.25;
    vertices[1] = 1; derivs[2] =  0.25; derivs[3] = -0.25;
    vertices[2] = 2; derivs[4] =  0.25; derivs[5] =  0.25;
    vertices[3] = 3; derivs[6] = -0.25; derivs[7] =  0.25;
    return;
  }
  
  const unsigned n4bit = 1<<0;
  const unsigned n5bit = 1<<1;
  const unsigned n6bit = 1<<2;
  const unsigned n7bit = 1<<3;
  
  vertices.resize(8);
  derivs.resize(16);
  msq_std::vector<size_t>::iterator v_iter = vertices.begin();
  msq_std::vector<double>::iterator d_iter = derivs.begin();
  
    // N_0
  if ((nodebits&(n4bit|n7bit)) != (n4bit|n7bit)) {  // if eiter adjacent mid-edge node is missing
    *v_iter = 0; ++v_iter;
    *d_iter = (nodebits&n7bit) ? 0.0 : -0.25; ++d_iter;
    *d_iter = (nodebits&n4bit) ? 0.0 : -0.25; ++d_iter;
  }
  
    // N_1
  if ((nodebits&(n4bit|n5bit)) != (n4bit|n5bit)) {  // if eiter adjacent mid-edge node is missing
    *v_iter = 1; ++v_iter;
    *d_iter = (nodebits&n5bit) ? 0.0 :  0.25; ++d_iter;
    *d_iter = (nodebits&n4bit) ? 0.0 : -0.25; ++d_iter;
  }
  
    // N_2
  if ((nodebits&(n5bit|n6bit)) != (n5bit|n6bit)) {  // if eiter adjacent mid-edge node is missing
    *v_iter = 2; ++v_iter;
    *d_iter = (nodebits&n5bit) ? 0.0 :  0.25; ++d_iter;
    *d_iter = (nodebits&n6bit) ? 0.0 :  0.25; ++d_iter;
  }
  
    // N_3
  if ((nodebits&(n6bit|n7bit)) != (n6bit|n7bit)) {  // if eiter adjacent mid-edge node is missing
    *v_iter = 3; ++v_iter;
    *d_iter = (nodebits&n7bit) ? 0.0 : -0.25; ++d_iter;
    *d_iter = (nodebits&n6bit) ? 0.0 :  0.25; ++d_iter;
  }
  
    // N_4
  if (nodebits&n4bit) {
    *v_iter = 4; ++v_iter;
    *d_iter =  0.0; ++d_iter;
    *d_iter = -0.5; ++d_iter;
  }
  
    // N_5
  if (nodebits&n5bit) {
    *v_iter = 5; ++v_iter;
    *d_iter =  0.5; ++d_iter;
    *d_iter =  0.0; ++d_iter;
  }
  
    // N_6
  if (nodebits&n6bit) {
    *v_iter = 6; ++v_iter;
    *d_iter =  0.0; ++d_iter;
    *d_iter =  0.5; ++d_iter;
  }
  
    // N_7
  if (nodebits&n7bit) {
    *v_iter = 7; ++v_iter;
    *d_iter = -0.5; ++d_iter;
    *d_iter =  0.0; ++d_iter;
  }
  
    // N_8 (mid-quad node) never contributes to Jacobian at element center!!!
  
    // Remove unused space
  vertices.resize( v_iter - vertices.begin() );
  derivs.resize( d_iter - derivs.begin() );  
}

} // namespace Mesquite
