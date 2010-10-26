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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Target2DUntangle.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DUntangleMu.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {


Target2DUntangleMu::~Target2DUntangleMu()
{}

std::string Target2DUntangleMu::get_name() const
  { return "untangle(" + mBaseMetric->get_name() + ")"; }

bool Target2DUntangleMu::evaluate( const MsqMatrix<2,2>& A, 
                                   const MsqMatrix<2,2>& W, 
                                   double& result, 
                                   MsqError& err )
{
  bool valid = mBaseMetric->evaluate( A, W, result, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
  const double d = mConstant - result;
  const double s = fabs(d) - d;
  result = 0.125*s*s*s;
  return true;
}

bool Target2DUntangleMu::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                             const MsqMatrix<2,2>& W,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_A,
                                             MsqError& err )
{
  bool valid = mBaseMetric->evaluate_with_grad( A, W, result, deriv_wrt_A, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
  if (mConstant < result) {
    const double s = result - mConstant;
    result = s*s*s;
    deriv_wrt_A *= 3*s*s;
  }
  else {
    result = 0;
    deriv_wrt_A = MsqMatrix<2,2>(0.0);
  }
  
  return true;
}


bool Target2DUntangleMu::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                             const MsqMatrix<2,2>& W,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_A,
                                             MsqMatrix<2,2> second_wrt_A[3],
                                             MsqError& err )
{
  bool valid = mBaseMetric->evaluate_with_hess( A, W, result, deriv_wrt_A, second_wrt_A, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
  if (mConstant < result) {
    const double s = result - mConstant;
    result = s*s*s;
    hess_scale( second_wrt_A, 3*s*s );
    pluseq_scaled_outer_product( second_wrt_A, 6*s, deriv_wrt_A );
    deriv_wrt_A *= 3*s*s;
  }
  else {
    result = 0;
    deriv_wrt_A = MsqMatrix<2,2>(0.0);
    second_wrt_A[0] = second_wrt_A[1] = second_wrt_A[2] = MsqMatrix<2,2>(0.0);
  }
  
  return true;
}

} // namespace MESQUITE_NS
