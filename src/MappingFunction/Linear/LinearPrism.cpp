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


/** \file LinearPrism.cpp
 *  \brief mapping function for linear prism
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "LinearPrism.hpp"

namespace Mesquite {

static const char* nonlinear_error 
 = "Attempt to use LinearPrism mapping function for a nonlinear element\n";


EntityTopology LinearPrism::element_topology() const
  { return PRISM; }
  
static const int edge_beg[] = { 0, 1, 2, 0, 1, 2, 3, 4, 5 };
static const int edge_end[] = { 1, 2, 0, 3, 4, 5, 4, 5, 3 };
static const int faces[5][5] = { { 4, 0, 1, 4, 3 },
                                 { 4, 1, 2, 5, 4 },
                                 { 4, 2, 0, 3, 5 },
                                 { 3, 0, 1, 2,-1 },
                                 { 3, 3, 4, 5,-1 } };

void LinearPrism::coefficients_at_corner( unsigned corner,
                                          unsigned nodebits,
                                          double* coeff_out,
                                          size_t& num_coeff,
                                          MsqError& err) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 6;
  coeff_out[0] = coeff_out[1] = coeff_out[2] = 0.0;
  coeff_out[3] = coeff_out[4] = coeff_out[5] = 0.0;
  coeff_out[corner] = 1.0;
}

void LinearPrism::coefficients_at_mid_edge( unsigned edge,
                                            unsigned nodebits,
                                            double* coeff_out,
                                            size_t& num_coeff,
                                            MsqError& err) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 6;
  coeff_out[0] = coeff_out[1] = coeff_out[2] = 0.0;
  coeff_out[3] = coeff_out[4] = coeff_out[5] = 0.0;
  coeff_out[edge_beg[edge]] = 0.5;
  coeff_out[edge_end[edge]] = 0.5;
}

void LinearPrism::coefficients_at_mid_face( unsigned face,
                                            unsigned nodebits,
                                            double* coeff_out,
                                            size_t& num_coeff,
                                            MsqError& err) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 6;
  coeff_out[0] = coeff_out[1] = coeff_out[2] = 0.0;
  coeff_out[3] = coeff_out[4] = coeff_out[5] = 0.0;
  const int n = faces[face][0];
  double f;
  if (n == 4) {
    f = 0.25;
    coeff_out[faces[face][4]] = f;
  }
  else
    f = MSQ_ONE_THIRD; 

  coeff_out[faces[face][1]] = f;
  coeff_out[faces[face][2]] = f;
  coeff_out[faces[face][3]] = f;
}

void LinearPrism::coefficients_at_mid_elem( unsigned nodebits,
                                            double* coeff_out,
                                            size_t& num_coeff,
                                            MsqError& err) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 6;
  const double sixth = 1.0/6.0;
  coeff_out[0] = sixth;
  coeff_out[1] = sixth;
  coeff_out[2] = sixth;
  coeff_out[3] = sixth;
  coeff_out[4] = sixth;
  coeff_out[5] = sixth;
}

void LinearPrism::derivatives_at_corner( unsigned corner, 
                                         unsigned nodebits,
                                         size_t* vertex_indices_out,
                                         double* d_coeff_d_xi_out,
                                         size_t& num_vtx,
                                         MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  int tri = (corner / 3); // 0 for xi=-1, 1 for xi=1
  int tv = corner % 3;    // offset of corner with xi=+/-1 triangle

  num_vtx = 4;
    // three vertices within the xi=+/-1 triangle
  vertex_indices_out[0] = 3*tri;
  vertex_indices_out[1] = 3*tri+1;
  vertex_indices_out[2] = 3*tri+2;
    // vertex adjacent to corner in other triangle
  vertex_indices_out[3] = 3 - 6*tri + corner;
  
    // three vertices within the xi=+/-1 triangle
  d_coeff_d_xi_out[ 0] =  0.0;
  d_coeff_d_xi_out[ 1] = -1.0;
  d_coeff_d_xi_out[ 2] = -1.0;
  d_coeff_d_xi_out[ 3] =  0.0;
  d_coeff_d_xi_out[ 4] =  1.0;
  d_coeff_d_xi_out[ 5] =  0.0;
  d_coeff_d_xi_out[ 6] =  0.0;
  d_coeff_d_xi_out[ 7] =  0.0;
  d_coeff_d_xi_out[ 8] =  1.0;
    // fix dxi value for input corner
  d_coeff_d_xi_out[3*tv] = tri - 0.5; 
    // vertex adjacent to corner in other triangle
  d_coeff_d_xi_out[ 9] =  0.5 - tri;
  d_coeff_d_xi_out[10] =  0.0;
  d_coeff_d_xi_out[11] =  0.0;
}


void LinearPrism::derivatives_at_mid_edge( unsigned edge, 
                                           unsigned nodebits,
                                           size_t* vertex_indices_out,
                                           double* d_coeff_d_xi_out,
                                           size_t& num_vtx,
                                           MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  
  int opp; // vertex opposite edge in same triagle
  
  switch (edge/3) {
    case 0:  // triangle at xi = -1
      opp = (edge+2)%3;

      num_vtx = 5;
        // vertices in this xi = -1 triagnle
      vertex_indices_out[0] = 0;
      vertex_indices_out[1] = 1;
      vertex_indices_out[2] = 2;
        // adjacent vertices in xi = 1 triangle
      vertex_indices_out[3] = 3 + edge;
      vertex_indices_out[4] = 3 + (edge+1)%3;

        // vertices in this xi = -1 triagnle
      d_coeff_d_xi_out[ 0] = -0.25;
      d_coeff_d_xi_out[ 1] = -1.00;
      d_coeff_d_xi_out[ 2] = -1.00;
      d_coeff_d_xi_out[ 3] = -0.25;
      d_coeff_d_xi_out[ 4] =  1.00;
      d_coeff_d_xi_out[ 5] =  0.00;
      d_coeff_d_xi_out[ 6] = -0.25;
      d_coeff_d_xi_out[ 7] =  0.00;
      d_coeff_d_xi_out[ 8] =  1.00;
        // clear dxi for vertex opposite edge in xi = -1 triangle
      d_coeff_d_xi_out[3*opp] = 0.0;
        // adjacent vertices in xi = 1 triangle
      d_coeff_d_xi_out[ 9] =  0.25;
      d_coeff_d_xi_out[10] =  0.00;
      d_coeff_d_xi_out[11] =  0.00;
      d_coeff_d_xi_out[12] =  0.25;
      d_coeff_d_xi_out[13] =  0.00;
      d_coeff_d_xi_out[14] =  0.00;
      break;

    case 1:  // lateral edges (not in either triangle)
      num_vtx = 6;
      vertex_indices_out[0] = 0;
      vertex_indices_out[1] = 1;
      vertex_indices_out[2] = 2;
      vertex_indices_out[3] = 3;
      vertex_indices_out[4] = 4;
      vertex_indices_out[5] = 5;
      
        // set all deta & dzeta values, zero all dxi values
      d_coeff_d_xi_out[ 0] =  0.0;
      d_coeff_d_xi_out[ 1] = -0.5;
      d_coeff_d_xi_out[ 2] = -0.5;
      d_coeff_d_xi_out[ 3] =  0.0;
      d_coeff_d_xi_out[ 4] =  0.5;
      d_coeff_d_xi_out[ 5] =  0.0;
      d_coeff_d_xi_out[ 6] =  0.0;
      d_coeff_d_xi_out[ 7] =  0.0;
      d_coeff_d_xi_out[ 8] =  0.5;
      d_coeff_d_xi_out[ 9] =  0.0;
      d_coeff_d_xi_out[10] = -0.5;
      d_coeff_d_xi_out[11] = -0.5;
      d_coeff_d_xi_out[12] =  0.0;
      d_coeff_d_xi_out[13] =  0.5;
      d_coeff_d_xi_out[14] =  0.0;
      d_coeff_d_xi_out[15] =  0.0;
      d_coeff_d_xi_out[16] =  0.0;
      d_coeff_d_xi_out[17] =  0.5;
      
        // set dxi values for end points of edge
      d_coeff_d_xi_out[3*(edge-3)] = -0.5;
      d_coeff_d_xi_out[3* edge   ] =  0.5;
      break;
    
    case 2:  // triangle at xi = 1
      opp = (edge+2)%3;

      num_vtx = 5;
        // vertices in this xi = 1 triagnle
      vertex_indices_out[0] = 3;
      vertex_indices_out[1] = 4;
      vertex_indices_out[2] = 5;
        // adjacent vertices in xi = 1 triangle
      vertex_indices_out[3] = edge - 6;
      vertex_indices_out[4] = (edge-5)%3;

        // vertices in this xi = -1 triagnle
      d_coeff_d_xi_out[ 0] =  0.25;
      d_coeff_d_xi_out[ 1] = -1.00;
      d_coeff_d_xi_out[ 2] = -1.00;
      d_coeff_d_xi_out[ 3] =  0.25;
      d_coeff_d_xi_out[ 4] =  1.00;
      d_coeff_d_xi_out[ 5] =  0.00;
      d_coeff_d_xi_out[ 6] =  0.25;
      d_coeff_d_xi_out[ 7] =  0.00;
      d_coeff_d_xi_out[ 8] =  1.00;
        // clear dxi for vertex opposite edge in xi = -1 triangle
      d_coeff_d_xi_out[3*opp] = 0.0;
        // adjacent vertices in xi = 1 triangle
      d_coeff_d_xi_out[ 9] = -0.25;
      d_coeff_d_xi_out[10] =  0.00;
      d_coeff_d_xi_out[11] =  0.00;
      d_coeff_d_xi_out[12] = -0.25;
      d_coeff_d_xi_out[13] =  0.00;
      d_coeff_d_xi_out[14] =  0.00;
      break;
  }
}
void LinearPrism::derivatives_at_mid_face( unsigned face, 
                                           unsigned nodebits,
                                           size_t* vertex_indices_out,
                                           double* d_coeff_d_xi_out,
                                           size_t& num_vtx,
                                           MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_vtx = 6;
  vertex_indices_out[0] = 0;
  vertex_indices_out[1] = 1;
  vertex_indices_out[2] = 2;
  vertex_indices_out[3] = 3;
  vertex_indices_out[4] = 4;
  vertex_indices_out[5] = 5;
  
    // dzeta is always zero for vertices 1 and 4
  d_coeff_d_xi_out[ 5] = 0.0;
  d_coeff_d_xi_out[14] = 0.0;
    // deta is always zero for vertices 2 and 5
  d_coeff_d_xi_out[ 7] = 0.0;
  d_coeff_d_xi_out[16] = 0.0;
  
  int opp; // start vtx of edge opposite from quad face
  int tri_offset; // offset in d_coeff_d_xi_out for triangle containing edge
  double sixth;
  
  if (face < 3) { // quad face
      // set all values
    d_coeff_d_xi_out[ 0] = -0.25;
    d_coeff_d_xi_out[ 1] = -0.50;
    d_coeff_d_xi_out[ 2] = -0.50;
    d_coeff_d_xi_out[ 3] = -0.25;
    d_coeff_d_xi_out[ 4] =  0.50;
    d_coeff_d_xi_out[ 5] =  0.00;
    d_coeff_d_xi_out[ 6] = -0.25;
    d_coeff_d_xi_out[ 7] =  0.00;
    d_coeff_d_xi_out[ 8] =  0.50;
    d_coeff_d_xi_out[ 9] =  0.25;
    d_coeff_d_xi_out[10] = -0.50;
    d_coeff_d_xi_out[11] = -0.50;
    d_coeff_d_xi_out[12] =  0.25;
    d_coeff_d_xi_out[13] =  0.50;
    d_coeff_d_xi_out[14] =  0.00;
    d_coeff_d_xi_out[15] =  0.25;
    d_coeff_d_xi_out[16] =  0.00;
    d_coeff_d_xi_out[17] =  0.50;
      // clear dxi for ends of edge opposite from face
    opp = (face+2)%3;
    d_coeff_d_xi_out[3*opp] = 0.0;
    d_coeff_d_xi_out[3*(opp+3)] = 0.0;
  }
  else { // triangular faces
      // set all xi values, zero all other values
    sixth = 1./6;
    d_coeff_d_xi_out[ 0] = -sixth;
    d_coeff_d_xi_out[ 1] =  0;
    d_coeff_d_xi_out[ 2] =  0;
    d_coeff_d_xi_out[ 3] = -sixth;
    d_coeff_d_xi_out[ 4] =  0;
    d_coeff_d_xi_out[ 5] =  0;
    d_coeff_d_xi_out[ 6] = -sixth;
    d_coeff_d_xi_out[ 7] =  0;
    d_coeff_d_xi_out[ 8] =  0;
    d_coeff_d_xi_out[ 9] =  sixth;
    d_coeff_d_xi_out[10] =  0;
    d_coeff_d_xi_out[11] =  0;
    d_coeff_d_xi_out[12] =  sixth;
    d_coeff_d_xi_out[13] =  0;
    d_coeff_d_xi_out[14] =  0;
    d_coeff_d_xi_out[15] =  sixth;
    d_coeff_d_xi_out[16] =  0;
    d_coeff_d_xi_out[17] =  0;
      // set deta and dzeta values for vertices in same triangle as edge
    tri_offset = 9 * (face - 3);  // either 0 or 9
    d_coeff_d_xi_out[tri_offset+1] = -1.0;
    d_coeff_d_xi_out[tri_offset+2] = -1.0;
    d_coeff_d_xi_out[tri_offset+4] =  1.0;
    d_coeff_d_xi_out[tri_offset+8] =  1.0;
  }
}
void LinearPrism::derivatives_at_mid_elem( unsigned nodebits,
                                           size_t* vertex_indices_out,
                                           double* d_coeff_d_xi_out,
                                           size_t& num_vtx,
                                           MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  const double sixth = 1./6;
  
  num_vtx = 6;;
  vertex_indices_out[0] = 0;
  vertex_indices_out[1] = 1;
  vertex_indices_out[2] = 2;
  vertex_indices_out[3] = 3;
  vertex_indices_out[4] = 4;
  vertex_indices_out[5] = 5;
  
  d_coeff_d_xi_out[ 0] = -sixth;
  d_coeff_d_xi_out[ 1] = -0.5;
  d_coeff_d_xi_out[ 2] = -0.5;
  d_coeff_d_xi_out[ 3] = -sixth;
  d_coeff_d_xi_out[ 4] =  0.5;
  d_coeff_d_xi_out[ 5] =  0.0;
  d_coeff_d_xi_out[ 6] = -sixth;
  d_coeff_d_xi_out[ 7] =  0.0;
  d_coeff_d_xi_out[ 8] =  0.5;
  d_coeff_d_xi_out[ 9] =  sixth;
  d_coeff_d_xi_out[10] = -0.5;
  d_coeff_d_xi_out[11] = -0.5;
  d_coeff_d_xi_out[12] =  sixth;
  d_coeff_d_xi_out[13] =  0.5;
  d_coeff_d_xi_out[14] =  0.0;
  d_coeff_d_xi_out[15] =  sixth;
  d_coeff_d_xi_out[16] =  0.0;
  d_coeff_d_xi_out[17] =  0.5;
}
  


} // namespace Mesquite
