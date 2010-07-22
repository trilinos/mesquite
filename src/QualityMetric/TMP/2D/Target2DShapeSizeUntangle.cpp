/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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
    (2010) jwfrank@sandia.gov   

  ***************************************************************** */


/** \file Target2DShapeSizeUntangle.cpp
 *  \brief 
 *  \author Jason Franks
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeUntangle.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {


Target2DShapeSizeUntangle::~Target2DShapeSizeUntangle()
{}

std::string Target2DShapeSizeUntangle::get_name() const
  { return "shape size untangle"; }

bool Target2DShapeSizeUntangle::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt(sqr_Frobenius(T) + 2.0*det(T));

  while (fabs(psi) < DBL_EPSILON) {
    T(0,0) += DBL_EPSILON;
    T(1,1) += DBL_EPSILON;
    frob_sqr = sqr_Frobenius(T);
    psi = sqrt( frob_sqr + 2.0*det(T) );
  }

  double mu = frob_sqr - 2.0*psi + 2.0;
  double d = (1 - mEps) - mu;
  double f = fabs(d) - d;
  result = f*f;
  return true;
}


bool Target2DShapeSizeUntangle::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                           const MsqMatrix<2,2>& W,
                                           double& result,
                                           MsqMatrix<2,2>& deriv_wrt_A,
                                           MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt(sqr_Frobenius(T) + 2.0*det(T));

  while (fabs(psi) < DBL_EPSILON) {
    T(0,0) += DBL_EPSILON;
    T(1,1) += DBL_EPSILON;
    frob_sqr = sqr_Frobenius(T);
    psi = sqrt( frob_sqr + 2.0*det(T) );
  }

  double mu = frob_sqr - 2.0*psi + 2.0;
  if (1 - mEps < mu) {
    double d = (1 - mEps) - mu;
    result = 4 * d*d;
    deriv_wrt_A = T;
    deriv_wrt_A *= (1.0 - 1.0/psi);
    deriv_wrt_A -= 1.0/psi * transpose_adj(T);
    deriv_wrt_A *= 2;
    deriv_wrt_A = 8 * d * (-1) * deriv_wrt_A;
    deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  }
  else {
    result = 0.0;
    deriv_wrt_A = MsqMatrix<2,2>(0.0);
  }
  return true;
}

} // namespace MESQUITE_NS
