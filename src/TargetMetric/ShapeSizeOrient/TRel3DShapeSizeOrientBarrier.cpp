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


/** \file TRel3DShapeSizeOrientBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel3DShapeSizeOrientBarrier.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

TRel3DShapeSizeOrientBarrier::~TRel3DShapeSizeOrientBarrier() {}

std::string TRel3DShapeSizeOrientBarrier::get_name() const
  { return "ShapeSizeOrientBarrier"; }

bool TRel3DShapeSizeOrientBarrier::evaluate( const MsqMatrix<3,3>& T, 
                                               double& result, 
                                               MsqError& err )
{
  double tau = det(T);
  if (invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  MsqMatrix<3,3> T_I(T);
  T_I(0,0) -= 1.0;
  T_I(1,1) -= 1.0;
  T_I(2,2) -= 1.0;
  result = sqr_Frobenius( T_I ) / (2*tau);
  return true;
}


bool TRel3DShapeSizeOrientBarrier::evaluate_with_grad( 
                                 const MsqMatrix<3,3>& T, 
                                 double& result, 
                                 MsqMatrix<3,3>& wrt_T,
                                 MsqError&  )
{
  double tau = det(T);
  if (invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  MsqMatrix<3,3> adjt = transpose_adj(T);
  MsqMatrix<3,3> T_I(T);
  T_I(0,0) -= 1.0;
  T_I(1,1) -= 1.0;
  T_I(2,2) -= 1.0;
  result = sqr_Frobenius(T_I) / (2*tau);

  wrt_T = T_I; // T - I
  wrt_T -= result * adjt;
  wrt_T *= 1.0/tau;
  return true;
}

bool TRel3DShapeSizeOrientBarrier::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                     double& result,
                                     MsqMatrix<3,3>& wrt_T,
                                     MsqMatrix<3,3> second[6],
                                     MsqError& err )
{
  double tau = det(T);
  if (invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  MsqMatrix<3,3> T_I(T);
  T_I(0,0) -= 1.0;
  T_I(1,1) -= 1.0;
  T_I(2,2) -= 1.0;
  result = sqr_Frobenius( T_I ) / (2*tau);

  MsqMatrix<3,3> adjt = transpose_adj(T);
  wrt_T = T_I; 
  wrt_T -= result * adjt;
  wrt_T *= 1.0/tau;
  
  set_scaled_outer_product( second, 2.0*result/(tau*tau), adjt );
  pluseq_scaled_sum_outer_product( second, -1.0/(tau*tau), T_I, adjt );
  pluseq_scaled_2nd_deriv_of_det( second, -result/tau, T );
  pluseq_scaled_I( second, 1.0/tau );
  
  return true;
}

} // namespace Mesquite
