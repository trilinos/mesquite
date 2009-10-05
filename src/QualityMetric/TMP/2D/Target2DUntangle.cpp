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


/** \file Target2DUntangle.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DUntangle.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {


Target2DUntangle::~Target2DUntangle()
{}

std::string Target2DUntangle::get_name() const
  { return "untangle"; }

bool Target2DUntangle::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  double tau = det(T);
  double d = tau - mGamma;
  double f = fabs(d) - d;
  result = f*f*f*f;
  return true;
}

bool Target2DUntangle::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                           const MsqMatrix<2,2>& W,
                                           double& result,
                                           MsqMatrix<2,2>& deriv_wrt_A,
                                           MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  double tau = det(T);
  if (tau < mGamma) {
    double d = tau - mGamma;
    result = 16 * d*d*d*d;
    deriv_wrt_A = 64 * d*d*d * transpose_adj(T);
    deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  }
  else {
    result = 0.0;
    deriv_wrt_A = MsqMatrix<2,2>(0.0);
  }
  return true;
}

bool Target2DUntangle::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                           const MsqMatrix<2,2>& W,
                                           double& result,
                                           MsqMatrix<2,2>& deriv_wrt_A,
                                           MsqMatrix<2,2> second_wrt_A[3],
                                           MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  double tau = det(T);
  if (tau < mGamma) {
    double d = tau - mGamma;
    result = 16 * d*d*d*d;
    const MsqMatrix<2,2> adjt = transpose_adj(T);
    deriv_wrt_A = 64 * d*d*d * adjt;
    deriv_wrt_A = deriv_wrt_A * transpose(Winv);
    set_scaled_outer_product( second_wrt_A, 192*d*d, adjt );
    pluseq_scaled_2nd_deriv_of_det( second_wrt_A, 64*d*d*d );
    second_deriv_wrt_product_factor( second_wrt_A, Winv );
  }
  else {
    result = 0.0;
    second_wrt_A[0] = second_wrt_A[1] = second_wrt_A[2] = deriv_wrt_A = MsqMatrix<2,2>(0.0);
  }
  return true;
}

} // namespace MESQUITE_NS
