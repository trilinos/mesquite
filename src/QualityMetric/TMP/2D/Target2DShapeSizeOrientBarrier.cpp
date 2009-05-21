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


/** \file Target2DShapeSizeOrientBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeOrientBarrier.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

msq_std::string Target2DShapeSizeOrientBarrier::get_name() const
  { return "ShapeSizeOrient"; }

bool Target2DShapeSizeOrientBarrier::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError& )
{
  MsqMatrix<2,2> T = A * inverse(W);
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  
  T(0,0) -= 1.0;
  T(1,1) -= 1.0;
  result = 0.5 * sqr_Frobenius(T) / d;
  return true;
}

bool Target2DShapeSizeOrientBarrier::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                               const MsqMatrix<2,2>& W,
                                               double& result,
                                               MsqMatrix<2,2>& deriv_wrt_A,
                                               MsqError& err )
{
  MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  
  MsqMatrix<2,2> D(T);
  D(0,0) -= 1.0;
  D(1,1) -= 1.0;
  double inv_d = 1.0/d;
  result = 0.5 * sqr_Frobenius(D) * inv_d;
  
  deriv_wrt_A = D;
  deriv_wrt_A -= result * transpose_adj(T);
  deriv_wrt_A *= inv_d;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  return true;
}

bool Target2DShapeSizeOrientBarrier::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                               const MsqMatrix<2,2>& W,
                                               double& result,
                                               MsqMatrix<2,2>& deriv_wrt_A,
                                               MsqMatrix<2,2> second_wrt_A[3],
                                               MsqError& err )
{
  MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  
  MsqMatrix<2,2> D(T);
  D(0,0) -= 1.0;
  D(1,1) -= 1.0;
  double inv_d = 1.0/d;
  result = 0.5 * sqr_Frobenius(D) * inv_d;
  
  MsqMatrix<2,2> adjt = transpose_adj(T);
  deriv_wrt_A = D;
  deriv_wrt_A -= result * adjt;
  deriv_wrt_A *= inv_d;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  set_scaled_outer_product( second_wrt_A, 2*result*inv_d*inv_d, adjt );
  pluseq_scaled_sum_outer_product( second_wrt_A, -inv_d*inv_d, D, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, -result * inv_d );
  pluseq_scaled_I( second_wrt_A, inv_d );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );
  
  return true;
}


} // namespace Mesquite
