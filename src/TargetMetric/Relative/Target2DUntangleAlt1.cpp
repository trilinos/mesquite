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


/** \file Target2DUntangleAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DUntangleAlt1.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {


Target2DUntangleAlt1::~Target2DUntangleAlt1()
{}

std::string Target2DUntangleAlt1::get_name() const
  { return "Untangle2"; }

bool Target2DUntangleAlt1::evaluate( const MsqMatrix<2,2>& A, 
                                     const MsqMatrix<2,2>& W, 
                                     double& result, 
                                     MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  double tau = det(T);
  result = 0.5 * (sqrt(tau*tau + mFactor) - tau);
  return true;
}

bool Target2DUntangleAlt1::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                               const MsqMatrix<2,2>& W,
                                               double& result,
                                               MsqMatrix<2,2>& deriv_wrt_A,
                                               MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  double tau = det(T);
  double g = sqrt(tau*tau + mFactor);
  double f = tau/g - 1;
  result = 0.5 * (g - tau);
  deriv_wrt_A = transpose_adj(T);
  deriv_wrt_A *= 0.5 * f;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  return true;
}

bool Target2DUntangleAlt1::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                               const MsqMatrix<2,2>& W,
                                               double& result,
                                               MsqMatrix<2,2>& deriv_wrt_A,
                                               MsqMatrix<2,2> second_wrt_A[3],
                                               MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  double tau = det(T);
  double g = sqrt(tau*tau + mFactor);
  double f = 0.5 * (tau/g - 1);
  result = 0.5 * (g - tau);
  
  deriv_wrt_A = adjt;
  deriv_wrt_A *= f;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  set_scaled_outer_product( second_wrt_A, 0.5*mFactor/(g*g*g), adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, f );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );
  
  return true;
}

} // namespace MESQUITE_NS
