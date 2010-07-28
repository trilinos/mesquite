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


/** \file Target3DShapeSizeUntangle.cpp
 *  \brief 
 *  \author Jason Franks
 */

#include "Mesquite.hpp"
#include "Target3DShapeSizeUntangle.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {


Target3DShapeSizeUntangle::~Target3DShapeSizeUntangle()
{}

std::string Target3DShapeSizeUntangle::get_name() const
  { return "shape size untangle"; }

bool Target3DShapeSizeUntangle::evaluate( const MsqMatrix<3,3>& A, 
                                          const MsqMatrix<3,3>& W, 
                                          double& result, 
                                          MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double nT = Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau-1;
  
  double mu = nT*nT*nT - 3*MSQ_SQRT_THREE*tau + mGamma*tau1*tau1;
  double d = (1 - mEps) - mu;
  double f = fabs(d) - d;
  result = f*f;
  return true;
}


bool Target3DShapeSizeUntangle::evaluate_with_grad( const MsqMatrix<3,3>& A,
                                                    const MsqMatrix<3,3>& W,
                                                    double& result,
                                                    MsqMatrix<3,3>& wrt_A,
                                                    MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  const double nT = Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau-1;
  double mu = nT*nT*nT - 3*MSQ_SQRT_THREE*tau + mGamma*tau1*tau1;

  if (1 - mEps < mu) {
    double d = (1 - mEps) - mu;
    result = 4 * d*d;
    wrt_A = T;
    wrt_A *= 3*nT;
    wrt_A -= (3*MSQ_SQRT_THREE - 2*mGamma*tau1) * transpose_adj(T);
    wrt_A = -8 * d * wrt_A;
    wrt_A = wrt_A * transpose(Winv);
  }
  else {
    result = 0.0;
    wrt_A = MsqMatrix<3,3>(0.0);
  }
  return true;
}

} // namespace MESQUITE_NS
