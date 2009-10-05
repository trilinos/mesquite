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


/** \file Target2DShapeBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeBarrier.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target2DShapeBarrier::get_name() const
  { return "ShapeBarrier"; }

bool Target2DShapeBarrier::evaluate( const MsqMatrix<2,2>& A, 
                                     const MsqMatrix<2,2>& W, 
                                     double& result, 
                                     MsqError& )
{
  const MsqMatrix<2,2> T = A * inverse(W);
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
    
  result = 0.5 * sqr_Frobenius(T) / d - 1;
  return true;
}

bool Target2DShapeBarrier::evaluate_with_grad( const MsqMatrix<2,2>& A,
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
  
  double inv_d = 1.0/d;
  result = 0.5 * sqr_Frobenius(T) * inv_d;
  deriv_wrt_A = T;
  deriv_wrt_A -= result * transpose_adj(T);
  deriv_wrt_A *= inv_d;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);

  result -= 1.0;
  return true;
}

/** \f$ \frac{\partial^2 \mu}{\partial T^2} 
      = \frac{1}{\tau} I_4 
      - \frac{1}{\tau^2} \left( T \otimes \frac{\partial \tau}{\partial T} 
                          + \frac{\partial \tau}{\partial T} \otimes T \right) 
      + \frac{|T|^2}{\tau^3} \left( \frac{\partial \tau}{\partial T} \otimes
                               \frac{\partial \tau}{\partial T} \right) 
      - \frac{|T|^2}{2 \tau^3} \frac{\partial^2 \tau}{\partial T^2} \f$
  */
bool Target2DShapeBarrier::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                               const MsqMatrix<2,2>& W,
                                               double& result,
                                               MsqMatrix<2,2>& deriv_wrt_A,
                                               MsqMatrix<2,2> second_wrt_A[3],
                                               MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
    
  double inv_d = 1.0/d;
  double f1 = sqr_Frobenius(T) * inv_d;
  result = 0.5 * f1;
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  deriv_wrt_A = T;
  deriv_wrt_A -= result * adjt;
  deriv_wrt_A *= inv_d;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  set_scaled_outer_product( second_wrt_A, f1 * inv_d * inv_d, adjt );
  pluseq_scaled_sum_outer_product( second_wrt_A, -inv_d*inv_d, T, adjt );
  pluseq_scaled_I( second_wrt_A, inv_d );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, -result * inv_d );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );

  result -= 1.0;
  return true;
}

} // namespace Mesquite
