/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
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


/** \file TargetMetric2D.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "TargetMetric2D.hpp"
#include "TargetMetricDimIndep.hpp"
#include "MsqError.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

TargetMetric2D::~TargetMetric2D() {}
     

bool TargetMetric2D::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                         const MsqMatrix<2,2>& W,
                                         double& result,
                                         MsqMatrix<2,2>& wrt_A,
                                         MsqError& err )
{
  bool valid = evaluate( A, W, result, err );
  if (MSQ_CHKERR(err) || !valid)
    return valid;
  
  wrt_A(0,0) = do_finite_difference( 0, 0, this, A, W, result, err ); MSQ_ERRZERO(err);
  wrt_A(0,1) = do_finite_difference( 0, 1, this, A, W, result, err ); MSQ_ERRZERO(err);
  wrt_A(1,0) = do_finite_difference( 1, 0, this, A, W, result, err ); MSQ_ERRZERO(err);
  wrt_A(1,1) = do_finite_difference( 1, 1, this, A, W, result, err ); MSQ_ERRZERO(err);
  return true;
}


bool TargetMetric2D::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                         const MsqMatrix<2,2>& W,
                                         double& result,
                                         MsqMatrix<2,2>& deriv_wrt_A,
                                         MsqMatrix<2,2> hess_wrt_A[3],
                                         MsqError& err )
{
  return do_numerical_hessian( this, A, W, result, deriv_wrt_A, hess_wrt_A, err );
}
  
} // namespace Mesquite
