/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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
 
    (2009) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file Target3DShapeSizeBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DShapeSizeBarrier.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target3DShapeSizeBarrier::get_name() const
  { return "ShapeSizeBarrier"; }

bool Target3DShapeSizeBarrier::evaluate( const MsqMatrix<3,3>& A, 
                                         const MsqMatrix<3,3>& W, 
                                         double& result, 
                                         MsqError&  )
{
  const MsqMatrix<3,3> T = A * inverse(W);
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const double nT = sqr_Frobenius(T);
  const double nadj = sqr_Frobenius(transpose_adj(T));
  const double f = 1/(tau*tau);
  result = nT + f*nadj - 6;
  return true;
}

bool Target3DShapeSizeBarrier::evaluate_with_grad( const MsqMatrix<3,3>& A,
                                                   const MsqMatrix<3,3>& W,
                                                   double& result,
                                                   MsqMatrix<3,3>& deriv_wrt_A,
                                                   MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double nadj = sqr_Frobenius(adjt);
  const double f = 1/(tau*tau);
  result = nT + f*nadj - 6;
  
  deriv_wrt_A = T;
  deriv_wrt_A *= (1+f*nT);
  deriv_wrt_A -= f * T * transpose(T) * T;
  deriv_wrt_A -= f/tau * nadj * adjt;
  deriv_wrt_A *= 2;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  return true;
}

bool Target3DShapeSizeBarrier::evaluate_with_hess( const MsqMatrix<3,3>& A,
                                                   const MsqMatrix<3,3>& W,
                                                   double& result,
                                                   MsqMatrix<3,3>& deriv_wrt_A,
                                                   MsqMatrix<3,3> second_wrt_A[6],
                                                   MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double nadj = sqr_Frobenius(adjt);
  const double f = 1/(tau*tau);
  result = nT + f*nadj - 6;
  
  //! \f$ \frac{\partial}{\partial T} |adj T|^2 \f$
  const MsqMatrix<3,3> dNadj_dT = 2 * (nT * T - T * transpose(T) * T);
  deriv_wrt_A = T;
  deriv_wrt_A -= f/tau * nadj * adjt;
  deriv_wrt_A *= 2;
  deriv_wrt_A += f * dNadj_dT;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);

    // calculate negative of 2nd wrt T of (|adj T|^2 / tau^2) (sec 3.2.2)
  set_scaled_2nd_deriv_norm_sqr_adj( second_wrt_A,    f,            T );
  pluseq_scaled_2nd_deriv_of_det(    second_wrt_A, -2*f*f*nadj*tau, T );
  pluseq_scaled_outer_product(       second_wrt_A,  6*f*f*nadj,     adjt );
  pluseq_scaled_sum_outer_product(   second_wrt_A, -2*f*f     *tau, adjt, dNadj_dT );
    // calculate 2nd wrt T of this metric
  pluseq_scaled_I( second_wrt_A, 2.0 );
    // calculate 2nd wrt A
  second_deriv_wrt_product_factor( second_wrt_A, Winv );

  return true;
}

} // namespace Mesquite
