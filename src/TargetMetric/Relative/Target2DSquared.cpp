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


/** \file Target2DSquared.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DSquared.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target2DSquared::get_name() const
  { return "sqr(" + mMetric->get_name() + ')'; }

bool Target2DSquared::evaluate( const MsqMatrix<2,2>& A, 
                                const MsqMatrix<2,2>& W, 
                                double& result, 
                                MsqError& err )
{
  bool rval = mMetric->evaluate( A, W, result, err );
  result *= result;
  return rval;
}

bool Target2DSquared::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                          const MsqMatrix<2,2>& W,
                                          double& result,
                                          MsqMatrix<2,2>& deriv_wrt_A,
                                          MsqError& err )
{
  bool rval = mMetric->evaluate_with_grad( A, W, result, deriv_wrt_A, err );
  deriv_wrt_A *= 2 * result;
  result *= result;
  return rval;
}

bool Target2DSquared::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                          const MsqMatrix<2,2>& W,
                                          double& result,
                                          MsqMatrix<2,2>& deriv_wrt_A,
                                          MsqMatrix<2,2> second_wrt_A[3],
                                          MsqError& err )
{
  bool rval = mMetric->evaluate_with_hess( A, W, result, deriv_wrt_A, second_wrt_A, err );
  second_wrt_A[0] *= 2 * result;
  second_wrt_A[1] *= 2 * result;
  second_wrt_A[2] *= 2 * result;
  pluseq_scaled_outer_product( second_wrt_A, 2.0, deriv_wrt_A );
  deriv_wrt_A *= 2 * result;
  result *= result;
  return rval;
}


} // namespace MESQUITE_NS
