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


/** \file TRel2DShapeSizeBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel2DShapeSizeBarrier.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string TRel2DShapeSizeBarrier::get_name() const
  { return "ShapeSizeBarrier"; }

bool TRel2DShapeSizeBarrier::evaluate( const MsqMatrix<2,2>& T, 
                                       double& result, 
                                       MsqError&  )
{
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

bool TRel2DShapeSizeBarrier::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                                 double& result,
                                                 MsqMatrix<2,2>& deriv_wrt_T,
                                                 MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
  
  deriv_wrt_T = T;
  deriv_wrt_T *= 2 + 2*f;
  deriv_wrt_T -= 2 * f/tau * nT * adjt;
  
  return true;
}

bool TRel2DShapeSizeBarrier::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                                 double& result,
                                                 MsqMatrix<2,2>& deriv_wrt_T,
                                                 MsqMatrix<2,2> second_wrt_T[3],
                                                 MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
  
  deriv_wrt_T = T;
  deriv_wrt_T *= 2 + 2*f;
  deriv_wrt_T -= 2 * f/tau * nT * adjt;

  set_scaled_sum_outer_product( second_wrt_T, -4*f/tau, T, adjt );
  pluseq_scaled_I( second_wrt_T, 2 + 2*f );
  pluseq_scaled_outer_product( second_wrt_T, 6*nT*f*f, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -2*nT*f/tau );

  return true;
}

} // namespace Mesquite
