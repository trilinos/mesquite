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


/** \file Target3DShapeOrientBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DShapeOrientBarrier.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target3DShapeOrientBarrier::get_name() const
  { return "ShapeOrientBarrier"; }

bool Target3DShapeOrientBarrier::evaluate( const MsqMatrix<3,3>& A, 
                                           const MsqMatrix<3,3>& W, 
                                           double& result, 
                                           MsqError&  )
{
  MsqMatrix<3,3> T = A * inverse(W);
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  result = 0.5/tau * (Frobenius( T ) - trace(T)/MSQ_SQRT_THREE);
  return true;
}

bool Target3DShapeOrientBarrier::evaluate_with_grad( const MsqMatrix<3,3>& A,
                                                     const MsqMatrix<3,3>& W,
                                                     double& result,
                                                     MsqMatrix<3,3>& deriv,
                                                     MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double norm = Frobenius(T);
  const double invroot = 1.0/MSQ_SQRT_THREE;
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  result = 0.5*inv_tau*(norm - invroot * trace(T));

  deriv = invnorm * T;
  deriv(0,0) -= invroot;
  deriv(1,1) -= invroot;
  deriv(2,2) -= invroot;
  deriv *= 0.5;
  deriv -= result * transpose_adj(T);
  deriv *= inv_tau;
  deriv = deriv * transpose(Winv);
  return true;
}


bool Target3DShapeOrientBarrier::evaluate_with_hess( const MsqMatrix<3,3>& A,
                                                     const MsqMatrix<3,3>& W,
                                                     double& result,
                                                     MsqMatrix<3,3>& deriv,
                                                     MsqMatrix<3,3> second[6],
                                                     MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double norm = Frobenius(T);
  const double invroot = 1.0/MSQ_SQRT_THREE;
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  const double f = norm - invroot * trace(T);
  result = 0.5 * inv_tau * f;

  const MsqMatrix<3,3> adjt = transpose_adj(T);
  deriv = invnorm * T;
  deriv(0,0) -= invroot;
  deriv(1,1) -= invroot;
  deriv(2,2) -= invroot;
  deriv *= 0.5;
  deriv -= result * adjt;
  deriv *= inv_tau;
  deriv = deriv * transpose(Winv);
  
  const double a = 0.5 * inv_tau * invnorm;
  set_scaled_outer_product( second, -a*invnorm*invnorm, T );
  pluseq_scaled_I( second, a );
  pluseq_scaled_outer_product( second, f*inv_tau*inv_tau*inv_tau, adjt );
  pluseq_scaled_2nd_deriv_of_det( second, -0.5*f*inv_tau*inv_tau, T );
  pluseq_scaled_sum_outer_product( second, -0.5*inv_tau*inv_tau*invnorm, T, adjt );
  pluseq_scaled_sum_outer_product_I( second, 0.5*inv_tau*inv_tau*invroot, adjt );
  second_deriv_wrt_product_factor( second, Winv );
  return true;
}


} // namespace Mesquite
