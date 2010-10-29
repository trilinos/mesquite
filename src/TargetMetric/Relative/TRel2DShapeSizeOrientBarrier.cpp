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


/** \file TRel2DShapeSizeOrientBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel2DShapeSizeOrientBarrier.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string TRel2DShapeSizeOrientBarrier::get_name() const
  { return "ShapeSizeOrientBarrier"; }

bool TRel2DShapeSizeOrientBarrier::evaluate( const MsqMatrix<2,2>& T, 
                                             double& result, 
                                             MsqError& )
{
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  
  MsqMatrix<2,2> T_I(T);
  T_I(0,0) -= 1.0;
  T_I(1,1) -= 1.0;
  result = 0.5 * sqr_Frobenius(T_I) / d;
  return true;
}

bool TRel2DShapeSizeOrientBarrier::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                                       double& result,
                                                       MsqMatrix<2,2>& deriv_wrt_T,
                                                       MsqError& err )
{
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
  
  deriv_wrt_T = D;
  deriv_wrt_T -= result * transpose_adj(T);
  deriv_wrt_T *= inv_d;
  
  return true;
}

bool TRel2DShapeSizeOrientBarrier::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                                       double& result,
                                                       MsqMatrix<2,2>& deriv_wrt_T,
                                                       MsqMatrix<2,2> second_wrt_T[3],
                                                       MsqError& err )
{
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
  deriv_wrt_T = D;
  deriv_wrt_T -= result * adjt;
  deriv_wrt_T *= inv_d;
  
  set_scaled_outer_product( second_wrt_T, 2*result*inv_d*inv_d, adjt );
  pluseq_scaled_sum_outer_product( second_wrt_T, -inv_d*inv_d, D, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -result * inv_d );
  pluseq_scaled_I( second_wrt_T, inv_d );
  
  return true;
}


} // namespace Mesquite
