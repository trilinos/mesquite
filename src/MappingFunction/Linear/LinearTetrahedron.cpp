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
#include "LinearTetrahedron.hpp"

namespace Mesquite {

static const char* nonlinear_error 
 = "Attempt to use LinearTetrahedron mapping function for a nonlinear element\n";
 
EntityTopology LinearTetrahedron::element_topology() const
  { return TETRAHEDRON; }

void LinearTetrahedron::coefficients_at_corner( unsigned corner,
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
  coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] = 0.0;
  coeff_out[corner] = 1.0;
}

void LinearTetrahedron::coefficients_at_mid_edge( unsigned edge,
                                                  unsigned nodebits,
                                                  double* coeff_out,
                                                  size_t& num_coeff,
                                                  MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  unsigned start_vtx, end_vtx, othr1_vtx, othr2_vtx;
  if (edge < 3) {
    start_vtx = edge;
    end_vtx   = (edge+1) % 3;
    othr1_vtx = (edge+2) % 3;
    othr2_vtx = 3;
  }
  else {
    start_vtx = edge-3;
    end_vtx = 3;
    othr1_vtx = (edge+1) % 3;
    othr2_vtx = (edge+2) % 3;
  }
  
  num_coeff = 4;
  coeff_out[start_vtx] = 0.5;
  coeff_out[  end_vtx] = 0.5;
  coeff_out[othr1_vtx] = 0.0;
  coeff_out[othr2_vtx] = 0.0;
}

void LinearTetrahedron::coefficients_at_mid_face( unsigned face,
                                                  unsigned nodebits,
                                                  double* coeff_out,
                                                  size_t& num_coeff,
                                                  MsqError& err ) const
{
  static const unsigned opposite_vtx[] = { 2, 0, 1, 3 };
  
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 4;
  coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] = MSQ_ONE_THIRD;
  coeff_out[opposite_vtx[face]] = 0.0;
}

void LinearTetrahedron::coefficients_at_mid_elem( unsigned nodebits,
                                                  double* coeff_out,
                                                  size_t& num_coeff,
                                                  MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  num_coeff = 4;
  coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] = 0.25;
}

static inline void tet_derivatives( unsigned nodebits,
                                    size_t* vertices,
                                    double* coeff_derivs,
                                    size_t& num_vtx,
                                    MsqError& err ) 
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
  }
  else {
    num_vtx = 4;
    vertices[0] = 0;
    vertices[1] = 1;
    vertices[2] = 2;
    vertices[3] = 3;
    
    coeff_derivs[ 0] = -1.0;
    coeff_derivs[ 1] = -1.0;
    coeff_derivs[ 2] = -1.0;

    coeff_derivs[ 3] =  1.0;
    coeff_derivs[ 4] =  0.0;
    coeff_derivs[ 5] =  0.0;

    coeff_derivs[ 6] =  0.0;
    coeff_derivs[ 7] =  1.0;
    coeff_derivs[ 8] =  0.0;

    coeff_derivs[ 9] =  0.0;
    coeff_derivs[10] =  0.0;
    coeff_derivs[11] =  1.0;
  }
}


void LinearTetrahedron::derivatives_at_corner( unsigned , 
                                               unsigned nodebits,
                                               size_t* vertex_indices_out,
                                               double* d_coeff_d_xi_out,
                                               size_t& num_vtx,
                                               MsqError& err ) const
{
  tet_derivatives( nodebits, vertex_indices_out, d_coeff_d_xi_out, num_vtx, err );
  MSQ_CHKERR(err);
}

void LinearTetrahedron::derivatives_at_mid_edge( unsigned , 
                                                 unsigned nodebits,
                                                 size_t* vertex_indices_out,
                                                 double* d_coeff_d_xi_out,
                                                 size_t& num_vtx,
                                                 MsqError& err ) const
{
  tet_derivatives( nodebits, vertex_indices_out, d_coeff_d_xi_out, num_vtx, err );
  MSQ_CHKERR(err);
}

void LinearTetrahedron::derivatives_at_mid_face( unsigned , 
                                                 unsigned nodebits,
                                                 size_t* vertex_indices_out,
                                                 double* d_coeff_d_xi_out,
                                                 size_t& num_vtx,
                                                 MsqError& err ) const
{
  tet_derivatives( nodebits, vertex_indices_out, d_coeff_d_xi_out, num_vtx, err );
  MSQ_CHKERR(err);
}

void LinearTetrahedron::derivatives_at_mid_elem( unsigned nodebits,
                                                 size_t* vertex_indices_out,
                                                 double* d_coeff_d_xi_out,
                                                 size_t& num_vtx,
                                                 MsqError& err ) const
{
  tet_derivatives( nodebits, vertex_indices_out, d_coeff_d_xi_out, num_vtx, err );
  MSQ_CHKERR(err);
}

} // namespace Mesquite
