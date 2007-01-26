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
                                          msq_std::vector<size_t>& indices,
                                          msq_std::vector<double>& derivs )
{
  indices.resize(5);
  indices[0] = 0;
  indices[1] = 1;
  indices[2] = 2;
  indices[3] = 3;
  indices[4] = 4;
    
  derivs.resize(15);
  
  derivs[ 0] = -value;
  derivs[ 1] = -value;
  derivs[ 2] = -0.125;
  
  derivs[ 3] =  value;
  derivs[ 4] = -value;
  derivs[ 5] = -0.125;
  
  derivs[ 6] =  value;
  derivs[ 7] =  value;
  derivs[ 8] = -0.125;
  
  derivs[ 9] = -value;
  derivs[10] =  value;
  derivs[11] = -0.125;
  
  derivs[12] =  0.0;
  derivs[13] =  0.0;
  derivs[14] =  0.5;
}

static inline void set_edge_derivatives( unsigned base_corner,
                                         double value,
                                         msq_std::vector<size_t>& indices,
                                         msq_std::vector<double>& derivs )
{
  const int direction = base_corner % 2;
  const int edge_beg =  base_corner;
  const int edge_end = (base_corner+1)%4;
  const int adj_end  = (base_corner+2)%4;
  const int adj_beg  = (base_corner+3)%4;
  const int dir_sign = 2*(edge_beg/2) - 1;
  const int oth_sign = 2*((edge_beg+1)/2%2) - 1;

  indices.resize(5);
  indices[0] = edge_beg;
  indices[1] = edge_end;
  indices[2] = adj_end;
  indices[3] = adj_beg;
  indices[4] = 4;

  derivs.resize(15);

  derivs[ 0+direction] =  2 * dir_sign * value;
  derivs[ 1-direction] =      oth_sign * value;
  derivs[ 2] = -0.25;

  derivs[ 3+direction] = -2 * dir_sign * value;
  derivs[ 4-direction] =      oth_sign * value;
  derivs[ 5] = -0.25;

  derivs[ 6+direction] =  0.0;
  derivs[ 7-direction] = -oth_sign * value;
  derivs[ 8]           =  0.0;

  derivs[ 9+direction] =  0.0;
  derivs[10-direction] = -oth_sign * value;
  derivs[11]           =  0.0;

  derivs[12] = 0.0;
  derivs[13] = 0.0;
  derivs[14] = 0.5;
}

static inline void set_corner_derivatives( unsigned corner,
                                           double value,
                                           msq_std::vector<size_t>& indices,
                                           msq_std::vector<double>& derivs )
{
  const unsigned adj_in_xi = (5 - corner) % 4;
  const unsigned adj_in_eta = 3 - corner;

  const int dxi_sign  = 2*((corner+1)/2%2)-1;
  const int deta_sign = 2*(corner/2) - 1;
  const double dxi_value = dxi_sign * value;
  const double deta_value = deta_sign * value;

  indices.resize(4);
  indices[0] = corner;
  indices[1] = adj_in_xi;
  indices[2] = adj_in_eta;
  indices[3] = 4;

  derivs.resize(12);

  derivs[ 0] =  dxi_value;
  derivs[ 1] =  deta_value;
  derivs[ 2] = -0.5;

  derivs[ 3] = -dxi_value;
  derivs[ 4] =  0.0;
  derivs[ 5] =  0.0;

  derivs[ 6] =  0.0;
  derivs[ 7] = -deta_value;
  derivs[ 8] =  0.0;

  derivs[ 9] =  0.0;
  derivs[10] =  0.0;
  derivs[11] =  0.5;
}

EntityTopology LinearPyramid::element_topology() const
  { return PYRAMID; }

void LinearPyramid::coefficients_at_corner( unsigned corner,
                                            unsigned nodebits,
                                            msq_std::vector<double>& coeff_out,
                                            MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  coeff_out.resize(5);
  coeff_out[0] = 0.0;
  coeff_out[1] = 0.0;
  coeff_out[2] = 0.0;
  coeff_out[3] = 0.0;
  coeff_out[4] = 0.0;
  coeff_out[corner] = 1.0;
}


void LinearPyramid::coefficients_at_mid_edge( unsigned edge,
                                            unsigned nodebits,
                                            msq_std::vector<double>& coeff_out,
                                            MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  coeff_out.resize(5);
  coeff_out[0] = 0.0;
  coeff_out[1] = 0.0;
  coeff_out[2] = 0.0;
  coeff_out[3] = 0.0;
  coeff_out[4] = 0.0;
  
  if (edge < 4) {
    coeff_out[edge] = 0.5;
    coeff_out[(edge+1)%4] = 0.5;
  }
  else {
    coeff_out[edge-4] = 0.5;
    coeff_out[4]      = 0.5;
  }
}

void LinearPyramid::coefficients_at_mid_face( unsigned face,
                                            unsigned nodebits,
                                            msq_std::vector<double>& coeff_out,
                                            MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }

  coeff_out.resize(5);
  if (face == 4) {
    coeff_out[0] = 0.25;
    coeff_out[1] = 0.25;
    coeff_out[2] = 0.25;
    coeff_out[3] = 0.25;
    coeff_out[4] = 0.00;
  }
  else {
    coeff_out[ face     ] = MSQ_ONE_THIRD;
    coeff_out[(face+1)%4] = MSQ_ONE_THIRD;
    coeff_out[(face+2)%4] = 0.0;
    coeff_out[(face+3)%4] = 0.0;
    coeff_out[      4   ] = MSQ_ONE_THIRD;
  }
}

void LinearPyramid::coefficients_at_mid_elem( 
                                            unsigned nodebits,
                                            msq_std::vector<double>& coeff_out,
                                            MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }

  coeff_out.resize(5);
  coeff_out[0] = 0.1875;
  coeff_out[1] = 0.1875;
  coeff_out[2] = 0.1875;
  coeff_out[3] = 0.1875;
  coeff_out[4] = 0.25;
}


void LinearPyramid::derivatives_at_corner( unsigned corner, 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  if (corner == 4) {
    set_equal_derivatives( 0.0, vertex_indices_out, d_coeff_d_xi_out );
  }
  else {
    set_corner_derivatives( corner, 0.5, vertex_indices_out, d_coeff_d_xi_out );
  }
}

void LinearPyramid::derivatives_at_mid_edge( unsigned edge, 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  if (edge < 4) {
    set_edge_derivatives( edge, 0.25, vertex_indices_out, d_coeff_d_xi_out );
  }
  else {
    set_corner_derivatives( edge-4, 0.25, vertex_indices_out, d_coeff_d_xi_out );
  }    
}


void LinearPyramid::derivatives_at_mid_face( unsigned face, 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  if (face == 4) {
    set_equal_derivatives( 0.25, vertex_indices_out, d_coeff_d_xi_out );
  }
  else {
    set_edge_derivatives( face, 1./6, vertex_indices_out, d_coeff_d_xi_out );
  }
}

void LinearPyramid::derivatives_at_mid_elem( 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  set_equal_derivatives( 0.1875, vertex_indices_out, d_coeff_d_xi_out );
}


} // namespace Mesquite
