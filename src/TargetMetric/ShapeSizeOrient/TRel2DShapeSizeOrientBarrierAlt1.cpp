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


/** \file TRel2DShapeSizeOrientBarrierAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel2DShapeSizeOrientBarrierAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string TRel2DShapeSizeOrientBarrierAlt1::get_name() const
  { return "ShapeSizeOrientBarrier1"; }

bool TRel2DShapeSizeOrientBarrierAlt1::evaluate( 
                                             const MsqMatrix<2,2>& T, 
                                             double& result, 
                                             MsqError& )
{
  double d = det(T);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  MsqMatrix<2,2> T_inv = 1/d * adj(T);
  T_inv(0,0) -= 1.0;
  T_inv(1,1) -= 1.0;
  result = sqr_Frobenius(T_inv);
  return true;
}

/** \f$ \frac{1}{\Tau^2}|T|^2 - \frac{2}{\Tau}tr(adj T) + 2 */
bool TRel2DShapeSizeOrientBarrierAlt1::evaluate_with_grad( 
                                             const MsqMatrix<2,2>& T,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_T,
                                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double it = 1.0/tau;
  result = it*(it*sqr_Frobenius(T) - 2.0*trace(T)) + 2.0;
  deriv_wrt_T = T;
  deriv_wrt_T *= it*it;
  deriv_wrt_T(0,0) -= it;
  deriv_wrt_T(1,1) -= it;
  deriv_wrt_T += it*it*(trace(T)-it*sqr_Frobenius(T))*adjt;
  deriv_wrt_T *= 2.0;
  return true;
}

bool TRel2DShapeSizeOrientBarrierAlt1::evaluate_with_hess( 
                                             const MsqMatrix<2,2>& T,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_T,
                                             MsqMatrix<2,2> second[3],
                                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double it = 1.0/tau;
  result = it*(it*sqr_Frobenius(T) - 2.0*trace(T)) + 2.0;
  deriv_wrt_T = T;
  deriv_wrt_T *= it*it;
  deriv_wrt_T(0,0) -= it;
  deriv_wrt_T(1,1) -= it;
  deriv_wrt_T += it*it*(trace(T)-it*sqr_Frobenius(T))*adjt;
  deriv_wrt_T *= 2.0;
  
  set_scaled_outer_product( second, it*it*it*(6*it*sqr_Frobenius(T) - 4*trace(T)), adjt );
  pluseq_scaled_I( second, 2*it*it );
  pluseq_scaled_2nd_deriv_of_det( second, 2*it*it*(trace(T) - it*sqr_Frobenius(T)) );
  pluseq_scaled_sum_outer_product( second, -4*it*it*it, T, adjt );
  pluseq_scaled_sum_outer_product_I( second, 2*it*it, adjt );
  
  return true;
}

} // namespace Mesquite
