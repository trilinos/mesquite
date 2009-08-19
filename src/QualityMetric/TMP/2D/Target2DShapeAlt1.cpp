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


/** \file Target2DShapeAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

msq_std::string Target2DShapeAlt1::get_name() const
  { return "Shape1"; }

bool Target2DShapeAlt1::evaluate( const MsqMatrix<2,2>& A, 
                              const MsqMatrix<2,2>& W, 
                              double& result, 
                              MsqError&  )
{
  MsqMatrix<2,2> T = A * inverse(W);
  const double tau = det(T);
  T = transpose(T) * T;
  T(0,0) -= tau;
  T(1,1) -= tau;
  result = sqr_Frobenius(T);
  return true;  
}

bool Target2DShapeAlt1::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                            const MsqMatrix<2,2>& W,
                                            double& result,
                                            MsqMatrix<2,2>& deriv_wrt_A,
                                            MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const MsqMatrix<2,2> TtT = transpose(T) * T;
  const double tau = det(T);
  const double nTtT = sqr_Frobenius(TtT);
  const double nT = sqr_Frobenius(T);
  result = nTtT + 2*tau*(tau - nT);
  
  deriv_wrt_A = T * TtT;
  deriv_wrt_A -= tau*T;
  deriv_wrt_A += (tau - 0.5*nT) * transpose_adj(T);
  deriv_wrt_A *= 4;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);

  return true;  
}

bool Target2DShapeAlt1::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                            const MsqMatrix<2,2>& W,
                                            double& result,
                                            MsqMatrix<2,2>& deriv_wrt_A,
                                            MsqMatrix<2,2> second_wrt_A[3],
                                            MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const MsqMatrix<2,2> TtT = transpose(T) * T;
  const double tau = det(T);
  const double nTtT = sqr_Frobenius(TtT);
  const double nT = sqr_Frobenius(T);
  result = nTtT + 2*tau*(tau - nT);
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  deriv_wrt_A = T * TtT;
  deriv_wrt_A -= tau*T;
  deriv_wrt_A += (tau - 0.5*nT) * adjt;
  deriv_wrt_A *= 4;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);

  set_scaled_outer_product( second_wrt_A, 1, T );
  second_wrt_A[1] = transpose(second_wrt_A[1]);
  second_wrt_A[0] += TtT;
  second_wrt_A[2] += TtT;
  const MsqMatrix<2,2> TTt = T * transpose(T);
  second_wrt_A[0](0,0) += TTt(0,0);
  second_wrt_A[0](1,1) += TTt(0,0);
  second_wrt_A[1](0,0) += TTt(0,1);
  second_wrt_A[1](1,1) += TTt(0,1);
  second_wrt_A[2](0,0) += TTt(1,1);
  second_wrt_A[2](1,1) += TTt(1,1);
  
  pluseq_scaled_I( second_wrt_A, -tau );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, -0.5*nT );
  pluseq_scaled_sum_outer_product( second_wrt_A, -1, T, adjt );
  
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, tau );
  pluseq_scaled_outer_product( second_wrt_A, 1, adjt );

  second_wrt_A[1] *= 4;
  second_wrt_A[0] *= 4;
  second_wrt_A[2] *= 4;
  
  second_deriv_wrt_product_factor( second_wrt_A, Winv );

  return true;
}

} // namespace Mesquite
