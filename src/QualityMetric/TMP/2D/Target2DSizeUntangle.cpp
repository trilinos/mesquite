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


/** \file Target2DSizeUntangle.cpp
 *  \brief 
 *  \author Jason Franks
 */

#include "Mesquite.hpp"
#include "Target2DSizeUntangle.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {


Target2DSizeUntangle::~Target2DSizeUntangle()
{}

std::string Target2DSizeUntangle::get_name() const
  { return "untangle"; }

bool Target2DSizeUntangle::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  double tau = det(T);
  double mu = (tau - 1.0)*(tau - 1.0);
  double d = (1.0 - mEps) - mu;
  double f = fabs(d) - d;
  result = f*f;
  return true;
}


bool  Target2DSizeUntangle::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                           const MsqMatrix<2,2>& W,
                                           double& result,
                                           MsqMatrix<2,2>& deriv_wrt_A,
                                           MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  double tau = det(T);
  double mu = (tau - 1.0)*(tau - 1.0);
  if (1 - mEps < mu) {
    double d = (1 - mEps) - mu;
    result = 4 * d*d;
    deriv_wrt_A = 16 * d * (1 - tau)*transpose_adj(T);
    deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  }
  else {
    result = 0.0;
    deriv_wrt_A = MsqMatrix<2,2>(0.0);
  }
  return true;
}

} // namespace MESQUITE_NS
