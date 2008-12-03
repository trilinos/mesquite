/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file InverseMeanRatio2D.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "InverseMeanRatio2D.hpp"
#include "MsqMatrix.hpp"

namespace Mesquite {

bool InverseMeanRatio2D::evaluate( const MsqMatrix<2,2>& A, 
                                   const MsqMatrix<2,2>& W, 
                                   double& result, 
                                   MsqError& err )
{
  const MsqMatrix<2,2> T = A * inverse(W);
  const double d = det( T );
  if (d < std::numeric_limits<double>::epsilon()) {
    result = 0.0;
    return false;
  }
  else {
    result = sqr_Frobenius(T) / (2 * d);\
    return true;
  }
}


bool InverseMeanRatio2D::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                             const MsqMatrix<2,2>& W,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_A,
                                             MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double d = det( T );
  if (d < std::numeric_limits<double>::epsilon()) {
    result = 0.0;
    deriv_wrt_A = MsqMatrix<2,2>(0.0);
    return false;
  }
  else {
    result = sqr_Frobenius(T) / (2 * d);
    deriv_wrt_A = transpose_adj(T);
    deriv_wrt_A *= -result;
    deriv_wrt_A += T;
    deriv_wrt_A *= 1.0/d;
    deriv_wrt_A = deriv_wrt_A * transpose(Winv);
    return true;
  }
}


bool InverseMeanRatio2D::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                             const MsqMatrix<2,2>& W,
                                             double& result,
                                             MsqMatrix<2,2>& dA,
                                             MsqMatrix<2,2> d2A[3],
                                             MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double d = det( T );
  if (d < std::numeric_limits<double>::epsilon()) {
    result = 0.0;
    dA = d2A[0] = d2A[1] = d2A[2] = MsqMatrix<2,2>(0.0);
    return false;
  }
  else {
    const double inv_det = 1.0/d;
    result = sqr_Frobenius(T) * 0.5 * inv_det;
    
    const MsqMatrix<2,2> AT = transpose_adj(T);
    dA = AT;
    dA *= -result;
    dA += T;
    dA *= inv_det;
    dA = dA * transpose(Winv);
    
    const double p3 = -result * inv_det;
    const double p1 = -2.0 * p3 * inv_det;
    const double p2 = -inv_det * inv_det;
    const MsqMatrix<2,2> AT_T_op_00 = outer( AT.row(0), T.row(0));
    const MsqMatrix<2,2> AT_T_op_11 = outer( AT.row(1), T.row(1));
    d2A[0] = p1 * outer( AT.row(0), AT.row(0))
           + p2 * (AT_T_op_00 + transpose(AT_T_op_00));
    d2A[1] = p1 * outer( AT.row(0), AT.row(1)) 
           + p2 * (outer( AT.row(0), T.row(1))
	     + outer( T.row(0), AT.row(1) ));
    d2A[2] = p1 * outer( AT.row(1), AT.row(1)) 
           + p2 * (AT_T_op_11 + transpose(AT_T_op_11));

    d2A[0](0,0) += inv_det;
    d2A[0](1,1) += inv_det;
    d2A[1](0,1) += p3;
    d2A[1](1,0) -= p3;
    d2A[2](0,0) += inv_det;
    d2A[2](1,1) += inv_det;
    
    d2A[0] = Winv * d2A[0] * transpose(Winv);
    d2A[1] = Winv * d2A[1] * transpose(Winv);
    d2A[2] = Winv * d2A[2] * transpose(Winv);
    
    return true;
  }
}



} // namespace Mesquite
