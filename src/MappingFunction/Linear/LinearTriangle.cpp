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

EntityTopology LinearTriangle::element_topology() const
  { return TRIANGLE; }

void LinearTriangle::coefficients( unsigned loc_dim,
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
  
  num_coeff = 3;
  switch (loc_dim) {
    case 0:
      coeff_out[0] = coeff_out[1] = coeff_out[2] = 0.0;
      coeff_out[loc_num] = 1.0;
      break;
    case 1:
      coeff_out[ loc_num     ] = 0.5;
      coeff_out[(loc_num+1)%3] = 0.5;
      coeff_out[(loc_num+2)%3] = 0.0;
      break;
    case 2:
      coeff_out[0] = coeff_out[1] = coeff_out[2] = MSQ_ONE_THIRD;
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

void LinearTriangle::derivatives( unsigned ,
                                  unsigned ,
                                  unsigned nodebits,
                                  size_t* vertices,
                                  MsqVector<2>* coeff_derivs,
                                  size_t& num_vtx,
                                  MsqError& err ) const
{
  if (nodebits) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
  }
  else {
    num_vtx = 3;
    vertices[0] = 0;
    vertices[1] = 1;
    vertices[2] = 2;
    
    coeff_derivs[0][0] = -1.0;
    coeff_derivs[0][1] = -1.0;
    coeff_derivs[1][0] = 1.0;
    coeff_derivs[1][1] = 0.0;
    coeff_derivs[2][0] = 0.0;
    coeff_derivs[2][1] = 1.0;
  }
}

} // namespace Mesquite
