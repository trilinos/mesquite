/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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


/** \file Target2DShapeSizeBarrierAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeBarrierAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target2DShapeSizeBarrierAlt1::get_name() const
  { return "ShapeSizeBarrierAlt1"; }

bool Target2DShapeSizeBarrierAlt1::evaluate( const MsqMatrix<2,2>& A, 
                                             const MsqMatrix<2,2>& W, 
                                             double& result, 
                                             MsqError&  )
{
  MsqMatrix<2,2> T = A * inverse(W);
  const double two_det = 2.0 * det(T);
  if (invalid_determinant(two_det)) { // barrier
    result = 0.0;
    return false;
  }
    
  const double frob_sqr = sqr_Frobenius(T);
  result = (frob_sqr - 2.0 * sqrt( frob_sqr + two_det ) + 2.0)/two_det;
  return true;
}

bool Target2DShapeSizeBarrierAlt1::evaluate_with_grad( const MsqMatrix<2,2>& A, 
                                                       const MsqMatrix<2,2>& W, 
                                                       double& result, 
                                                       MsqMatrix<2,2>& deriv_wrt_A,
                                                       MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  const double frob_sqr = sqr_Frobenius(T);
  const double psi = sqrt( frob_sqr + 2.0*det(T) );
  const double v = frob_sqr - 2.0 * psi + 2.0;
  result = v / (2*d);

    // deriv of V wrt T
  MsqMatrix<2,2> adjt = transpose_adj(T);
  MsqMatrix<2,2> v_wrt_T(T);
  v_wrt_T *= (1.0 - 1.0/psi);
  v_wrt_T -= 1.0/psi * adjt;
  v_wrt_T *= 2;
  
    // deriv of mu wrt T
  deriv_wrt_A = v_wrt_T;
  deriv_wrt_A *= 0.5/d;
  deriv_wrt_A -= v / (2*d*d) * adjt;
  
    // deriv of mu wrt A
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  return true;
}

bool Target2DShapeSizeBarrierAlt1::evaluate_with_hess( const MsqMatrix<2,2>& A, 
                                                       const MsqMatrix<2,2>& W, 
                                                       double& result, 
                                                       MsqMatrix<2,2>& deriv_wrt_A,
                                                       MsqMatrix<2,2> second[3],
                                                       MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  const double frob_sqr = sqr_Frobenius(T);
  const double psi = sqrt( frob_sqr + 2.0*det(T) );
  const double v = frob_sqr - 2.0 * psi + 2.0;
  result = v / (2*d);

    // deriv of V wrt T
  MsqMatrix<2,2> adjt = transpose_adj(T);
  MsqMatrix<2,2> v_wrt_T(T);
  v_wrt_T *= (1.0 - 1.0/psi);
  v_wrt_T -= 1.0/psi * adjt;
  v_wrt_T *= 2;
  
    // deriv of mu wrt T
  deriv_wrt_A = v_wrt_T;
  deriv_wrt_A *= 0.5/d;
  deriv_wrt_A -= v / (2*d*d) * adjt;
  
    // deriv of mu wrt A
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
    // second of V wrt T 
  const double s = T(0,1) - T(1,0);
  const double tr = trace(T);
  const double f = -2.0/(psi*psi*psi);
  second[0](0,0) = second[1](0,1) = second[2](1,1) =  f*s*s;
  second[0](0,1) = second[0](1,0) = second[1](1,1) = -f*s*tr;
  second[1](0,0) = second[2](0,1) = second[2](1,0) =  f*s*tr;
  second[0](1,1) = second[2](0,0) = -(second[1](1,0) = -f*tr*tr);
  pluseq_scaled_I( second, 2 );
  
    // second of mu wrt T 
  const double x = 1.0/(2*d);
  second[0] *= x;
  second[1] *= x;
  second[2] *= x;
  pluseq_scaled_2nd_deriv_of_det( second, v/(-2*d*d) );
  pluseq_scaled_outer_product( second, v/(d*d*d), adjt );
  pluseq_scaled_sum_outer_product( second, -1/(2*d*d), v_wrt_T, adjt );
  
    // second of mu wrt A
  second_deriv_wrt_product_factor( second, Winv );
  
  return true;
}

} // namespace Mesquite
