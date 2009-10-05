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
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target2DShapeSizeBarrierAlt1::get_name() const
  { return "ShapeSizeBarrier"; }

bool Target2DShapeSizeBarrierAlt1::evaluate( const MsqMatrix<2,2>& A, 
                                             const MsqMatrix<2,2>& W, 
                                             double& result, 
                                             MsqError&  )
{
  const MsqMatrix<2,2> T = A * inverse(W);
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
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
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
  
  deriv_wrt_A = T;
  deriv_wrt_A *= 2 + 2*f;
  deriv_wrt_A -= 2 * f/tau * nT * adjt;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  return true;
}

bool Target2DShapeSizeBarrierAlt1::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                                       const MsqMatrix<2,2>& W,
                                                       double& result,
                                                       MsqMatrix<2,2>& deriv_wrt_A,
                                                       MsqMatrix<2,2> second_wrt_A[3],
                                                       MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
  
  deriv_wrt_A = T;
  deriv_wrt_A *= 2 + 2*f;
  deriv_wrt_A -= 2 * f/tau * nT * adjt;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
/*  
    // calculate negative of 2nd wrt T of (|adj T|^2 / tau^2) (eqn 3.75)
  set_scaled_2nd_deriv_norm_sqr_adj( second_wrt_A, -f, T );
  pluseq_scaled_outer_product( second_wrt_A, -6 * f*f * nT, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, 2*f*f*tau*nT );
  pluseq_scaled_sum_outer_product( second_wrt_A, 4*f*f*tau, adjt, T );
    // calculate 2nd wrt T of this metric
  pluseq_scaled_I( second_wrt_A, 2.0 );
    // calculate 2nd wrt A
  second_deriv_wrt_product_factor( second_wrt_A, Winv );
*/
  set_scaled_sum_outer_product( second_wrt_A, -4*f/tau, T, adjt );
  pluseq_scaled_I( second_wrt_A, 2 + 2*f );
  pluseq_scaled_outer_product( second_wrt_A, 6*nT*f*f, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, -2*nT*f/tau );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );

  return true;
}

} // namespace Mesquite
