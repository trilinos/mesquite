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

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "LinearTriangle.hpp"

namespace Mesquite {

static const char* nonlinear_error 
 = "Attempt to use LinearTriangle mapping function for a nonlinear element\n";
 
static const char* dimension_error
 = "Cannot do midface evaluation for 2D elements.\n";

EntityTopology LinearTriangle::element_topology() const
  { return TRIANGLE; }

void LinearTriangle::coefficients_at_corner( unsigned corner,
                                             unsigned nodebits,
                                             msq_std::vector<double>& coeff_out,
                                             MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  coeff_out.resize(3);
  coeff_out[0] = coeff_out[1] = coeff_out[2] = 0.0;
  coeff_out[corner] = 1.0;
}

void LinearTriangle::coefficients_at_mid_edge( unsigned edge,
                                               unsigned nodebits,
                                               msq_std::vector<double>& coeff_out,
                                               MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  coeff_out.resize(3);
  const unsigned start_vtx =  edge;
  const unsigned   end_vtx = (edge+1)%3;
  const unsigned other_vtx = (edge+2)%3;
  coeff_out[start_vtx] = 0.5;
  coeff_out[  end_vtx] = 0.5;
  coeff_out[other_vtx] = 0.0;
}

void LinearTriangle::coefficients_at_mid_face( unsigned ,
                                               unsigned ,
                                               msq_std::vector<double>& ,
                                               MsqError& err ) const
{
  MSQ_SETERR(err)(dimension_error, MsqError::UNSUPPORTED_ELEMENT );
  return;
}

void LinearTriangle::coefficients_at_mid_elem( unsigned nodebits,
                                               msq_std::vector<double>& coeff_out,
                                               MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  coeff_out.resize(3);
  coeff_out[0] = coeff_out[1] = coeff_out[2] = MSQ_ONE_THIRD;
}

static inline void triangle_derivatives( unsigned nodebits,
                                         msq_std::vector<size_t>& vertices,
                                         msq_std::vector<double>& coeff_derivs,
                                         MsqError& err ) 
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
  }
  else {
    vertices.resize(3);
    vertices[0] = 0;
    vertices[1] = 1;
    vertices[2] = 2;
    
    coeff_derivs.resize(6);
    coeff_derivs[0] = -1.0;
    coeff_derivs[1] = -1.0;
    coeff_derivs[2] = 1.0;
    coeff_derivs[3] = 0.0;
    coeff_derivs[4] = 0.0;
    coeff_derivs[5] = 1.0;
  }
}


void LinearTriangle::derivatives_at_corner( unsigned , 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const
{
  triangle_derivatives( nodebits, vertex_indices_out, d_coeff_d_xi_out, err );
  MSQ_CHKERR(err);
}

void LinearTriangle::derivatives_at_mid_edge( unsigned , 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const
{
  triangle_derivatives( nodebits, vertex_indices_out, d_coeff_d_xi_out, err );
  MSQ_CHKERR(err);
}

void LinearTriangle::derivatives_at_mid_face( unsigned , 
                              unsigned ,
                              msq_std::vector<size_t>& ,
                              msq_std::vector<double>& ,
                              MsqError& err ) const
{
  MSQ_SETERR(err)(dimension_error, MsqError::UNSUPPORTED_ELEMENT );
  return;
}

void LinearTriangle::derivatives_at_mid_elem(
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const
{
  triangle_derivatives( nodebits, vertex_indices_out, d_coeff_d_xi_out, err );
  MSQ_CHKERR(err);
}

} // namespace Mesquite
