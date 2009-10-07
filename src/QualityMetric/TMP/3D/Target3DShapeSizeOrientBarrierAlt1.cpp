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


/** \file Target3DShapeSizeOrientBarrierAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DShapeSizeOrientBarrierAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target3DShapeSizeOrientBarrierAlt1::get_name() const
  { return "ShapeSizeOrientBarrier1"; }

bool Target3DShapeSizeOrientBarrierAlt1::evaluate( 
                                             const MsqMatrix<3,3>& A, 
                                             const MsqMatrix<3,3>& W, 
                                             double& result, 
                                             MsqError& )
{
  double d = det(A);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  MsqMatrix<3,3> T_inv = 1/d * W * adj(A);
  T_inv(0,0) -= 1.0;
  T_inv(1,1) -= 1.0;
  T_inv(2,2) -= 1.0;
  result = sqr_Frobenius(T_inv);
  return true;
}

/** \f$ \frac{1}{\Tau^2}|adj T|^2 - \frac{2}{\Tau}tr(adj T) + 3 */
bool Target3DShapeSizeOrientBarrierAlt1::evaluate_with_grad( 
                                             const MsqMatrix<3,3>& A,
                                             const MsqMatrix<3,3>& W,
                                             double& result,
                                             MsqMatrix<3,3>& deriv_wrt_A,
                                             MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<3,3> adjt = adj(T);
  const double it = 1.0/tau;
  result = it*(it*sqr_Frobenius(adjt) - 2.0*trace(adjt)) + 3.0;
  
  deriv_wrt_A = T;
  deriv_wrt_A *= sqr_Frobenius(T);
  deriv_wrt_A -= T * transpose(T) * T;
  deriv_wrt_A *= it*it;
  
  deriv_wrt_A += it*it*(trace(adjt)-it*sqr_Frobenius(adjt))*transpose(adjt);

  double f = trace(T) * it;
  deriv_wrt_A(0,0) -= f;
  deriv_wrt_A(1,1) -= f;
  deriv_wrt_A(2,2) -= f;
  
  deriv_wrt_A += it*transpose(T);

  deriv_wrt_A *= 2.0;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  return true;
}

bool Target3DShapeSizeOrientBarrierAlt1::evaluate_with_hess( 
                                             const MsqMatrix<3,3>& A,
                                             const MsqMatrix<3,3>& W,
                                             double& result,
                                             MsqMatrix<3,3>& deriv_wrt_A,
                                             MsqMatrix<3,3> second[6],
                                             MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<3,3> adjt = adj(T);
  const double it = 1.0/tau;
  const double nadjt = sqr_Frobenius(adjt);
  const double nT = sqr_Frobenius(T);
  const double tadjT = trace(adjt);
  result = it*(it*nadjt - 2.0*tadjT) + 3.0;
  
  const MsqMatrix<3,3> TTtT = T * transpose(T) * T;
  deriv_wrt_A = T;
  deriv_wrt_A *= nT;
  deriv_wrt_A -= TTtT;
  deriv_wrt_A *= it*it;
 
  deriv_wrt_A += it*it*(tadjT-it*nadjt)*transpose(adjt);

  const double tT = trace(T);
  double f = tT * it;
  deriv_wrt_A(0,0) -= f;
  deriv_wrt_A(1,1) -= f;
  deriv_wrt_A(2,2) -= f;
  
  deriv_wrt_A += it*transpose(T);

  deriv_wrt_A *= 2.0;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  set_scaled_2nd_deriv_norm_sqr_adj( second, it*it, T );
/*  
  const double yf = -it*it*it*it;
  pluseq_scaled_2nd_deriv_of_det( second, yf*2*nadjt*tau, T );
  pluseq_scaled_outer_product( second, yf*2*nadjt, transpose(adjt) );
  pluseq_scaled_outer_product( second, yf*-8*nadjt, transpose(adjt) );
  MsqMatrix<3,3> dnadj_dT = 2 * (sqr_Frobenius(T) * T - T * transpose(T) * T);
  pluseq_scaled_sum_outer_product( second, yf * 2 * tau, dnadj_dT, transpose(adjt) );
  const double sf = -2;
  pluseq_scaled_2nd_deriv_tr_adj( second, sf * it );
  const double zf = -it*it*sf;
  pluseq_scaled_2nd_deriv_of_det( second, zf * trace(adjt), T );
  pluseq_scaled_outer_product( second, zf * -2 * trace(adjt) * it, transpose(adjt) );
  MsqMatrix<3,3> dtradj_dT = -transpose(T);
  dtradj_dT(0,0) += trace(T);
  dtradj_dT(1,1) += trace(T);
  dtradj_dT(2,2) += trace(T);
  pluseq_scaled_sum_outer_product( second, zf, dtradj_dT, transpose(adjt) );
*/
  const double yf = -it*it*it*it;
  const double sf = -2;
  const double zf = -it*it*sf;

  pluseq_scaled_2nd_deriv_of_det( second, yf*2*nadjt*tau + zf*tadjT, T );
  pluseq_scaled_outer_product( second, yf*-6*nadjt - 2*zf*tadjT*it, transpose(adjt) );
  MsqMatrix<3,3> dnadj_dT = 2 * (nT * T - TTtT);
  pluseq_scaled_sum_outer_product( second, yf * 2 * tau, dnadj_dT, transpose(adjt) );
  pluseq_scaled_2nd_deriv_tr_adj( second, sf * it );
  MsqMatrix<3,3> dtradj_dT = -transpose(T);
  dtradj_dT(0,0) += tT;
  dtradj_dT(1,1) += tT;
  dtradj_dT(2,2) += tT;
  pluseq_scaled_sum_outer_product( second, zf, dtradj_dT, transpose(adjt) );

  second_deriv_wrt_product_factor( second, Winv );

  return true;
}

} // namespace Mesquite
