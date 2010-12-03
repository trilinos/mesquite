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


/** \file TRel2DUntangle.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel2DUntangleMu.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {


TRel2DUntangleMu::~TRel2DUntangleMu()
{}

std::string TRel2DUntangleMu::get_name() const
  { return "untangle(" + mBaseMetric->get_name() + ")"; }

bool TRel2DUntangleMu::evaluate( const MsqMatrix<2,2>& T, 
                                 double& result, 
                                 MsqError& err )
{
  bool valid = mBaseMetric->evaluate( T, result, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
  const double d = mConstant - result;
  const double s = fabs(d) - d;
  result = 0.125*s*s*s;
  return true;
}

bool TRel2DUntangleMu::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                           double& result,
                                           MsqMatrix<2,2>& deriv_wrt_T,
                                           MsqError& err )
{
  bool valid = mBaseMetric->evaluate_with_grad( T, result, deriv_wrt_T, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
  if (mConstant < result) {
    const double s = result - mConstant;
    result = s*s*s;
    deriv_wrt_T *= 3*s*s;
  }
  else {
    result = 0;
    deriv_wrt_T = MsqMatrix<2,2>(0.0);
  }
  
  return true;
}


bool TRel2DUntangleMu::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                           double& result,
                                           MsqMatrix<2,2>& deriv_wrt_T,
                                           MsqMatrix<2,2> second_wrt_T[3],
                                           MsqError& err )
{
  bool valid = mBaseMetric->evaluate_with_hess( T, result, deriv_wrt_T, second_wrt_T, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
  if (mConstant < result) {
    const double s = result - mConstant;
    result = s*s*s;
    hess_scale( second_wrt_T, 3*s*s );
    pluseq_scaled_outer_product( second_wrt_T, 6*s, deriv_wrt_T );
    deriv_wrt_T *= 3*s*s;
  }
  else {
    result = 0;
    deriv_wrt_T = MsqMatrix<2,2>(0.0);
    second_wrt_T[0] = second_wrt_T[1] = second_wrt_T[2] = MsqMatrix<2,2>(0.0);
  }
  
  return true;
}

} // namespace MESQUITE_NS
