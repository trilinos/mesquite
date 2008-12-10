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

void LinearQuadrilateral::coefficients_at_corner( unsigned corner,
                                                  unsigned nodebits,
                                                  double* coeff_out,
                                                  size_t& num_coeff,
                                                  MsqError& err) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 4;
  coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] = 0.0;
  coeff_out[corner] = 1.0;
}

void LinearQuadrilateral::coefficients_at_mid_edge( unsigned edge,
                                                    unsigned nodebits,
                                                    double* coeff_out,
                                                    size_t& num_coeff,
                                                    MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  const unsigned start_vtx = edge;
  const unsigned   end_vtx = (edge+1)%4;
  const unsigned othr1_vtx = (edge+2)%4;
  const unsigned othr2_vtx = (edge+3)%4;
  
  num_coeff = 4;
  coeff_out[start_vtx] = coeff_out[  end_vtx] = 0.5;
  coeff_out[othr1_vtx] = coeff_out[othr2_vtx] = 0.0;
}

void LinearQuadrilateral::coefficients_at_mid_face( unsigned ,
                                                    unsigned ,
                                                    double* ,
                                                    size_t& ,
                                                    MsqError& err) const
{
  MSQ_SETERR(err)("Request for mid-face mapping function value"
                  "for a quadrilateral element", MsqError::UNSUPPORTED_ELEMENT);
}

void LinearQuadrilateral::coefficients_at_mid_elem( unsigned nodebits,
                                                    double* coeff_out,
                                                    size_t& num_coeff,
                                                    MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 4;
  coeff_out[0] = 0.25;
  coeff_out[1] = 0.25;
  coeff_out[2] = 0.25;
  coeff_out[3] = 0.25;
}

const unsigned  xi = 0;
const unsigned eta = 1;
const int sign[2][4] = {{ -1,  1,  1, -1 },  // xi
                        { -1, -1,  1,  1 }}; // eta

void LinearQuadrilateral::derivatives_at_corner( unsigned corner, 
                                                 unsigned nodebits,
                                                 size_t* vertex_indices,
                                                 double* derivs,
                                                 size_t& num_vtx,
                                                 MsqError& err  ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  const unsigned adj_in_xi = (5 - corner) % 4;
  const unsigned adj_in_eta = 3 - corner;
  
  num_vtx = 3;
  vertex_indices[0] = corner;
  vertex_indices[1] = adj_in_xi;
  vertex_indices[2] = adj_in_eta;
  
  derivs[0] = 0.5 * sign[ xi][corner];
  derivs[1] = 0.5 * sign[eta][corner];
  derivs[2] = 0.5 * sign[ xi][adj_in_xi ];
  derivs[3] = 0.0;
  derivs[4] = 0.0;
  derivs[5] = 0.5 * sign[eta][adj_in_eta];
}

void LinearQuadrilateral::derivatives_at_mid_edge( unsigned edge, 
                                                   unsigned nodebits,
                                                   size_t* vertices,
                                                   double* derivs,
                                                   size_t& num_vtx,
                                                   MsqError&  err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
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
  
  derivs[2*start_vtx + direction] = 0.5 * sign[direction][start_vtx];
  derivs[2*  end_vtx + direction] = 0.5 * sign[direction][  end_vtx];
  derivs[2*othr1_vtx + direction] = 0.0;
  derivs[2*othr2_vtx + direction] = 0.0;
 
  derivs[0 + orthogonal] = 0.25 * sign[orthogonal][0];
  derivs[2 + orthogonal] = 0.25 * sign[orthogonal][1];
  derivs[4 + orthogonal] = 0.25 * sign[orthogonal][2];
  derivs[6 + orthogonal] = 0.25 * sign[orthogonal][3];
}

  
void LinearQuadrilateral::derivatives_at_mid_face( unsigned , 
                                                   unsigned ,
                                                   size_t* ,
                                                   double* ,
                                                   size_t& ,
                                                   MsqError& err ) const

{
  MSQ_SETERR(err)("Request for mid-face mapping function derivative"
                  "for a quadrilateral element", MsqError::UNSUPPORTED_ELEMENT);
}


void LinearQuadrilateral::derivatives_at_mid_elem( unsigned nodebits,
                                                   size_t* vertices,
                                                   double* derivs,
                                                   size_t& num_vtx,
                                                   MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_vtx = 4;
  vertices[0] = 0; derivs[0] = -0.25; derivs[1] = -0.25;
  vertices[1] = 1; derivs[2] =  0.25; derivs[3] = -0.25;
  vertices[2] = 2; derivs[4] =  0.25; derivs[5] =  0.25;
  vertices[3] = 3; derivs[6] = -0.25; derivs[7] =  0.25;
}

} // namespace Mesquite
