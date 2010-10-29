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


/** \file TRel3DMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "TRel3DMetric.hpp"
#include "TargetMetricDimIndep.hpp"

namespace MESQUITE_NS {

TRel3DMetric::~TRel3DMetric() {}
     

bool TRel3DMetric::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                       double& result,
                                       MsqMatrix<3,3>& wrt_T,
                                       MsqError& err )
{
  bool valid = evaluate( T, result, err );
  if (MSQ_CHKERR(err) || !valid)
    return valid;
  
  wrt_T(0,0) = do_finite_difference( 0, 0, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(0,1) = do_finite_difference( 0, 1, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(0,2) = do_finite_difference( 0, 2, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(1,0) = do_finite_difference( 1, 0, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(1,1) = do_finite_difference( 1, 1, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(1,2) = do_finite_difference( 1, 2, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(2,0) = do_finite_difference( 2, 0, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(2,1) = do_finite_difference( 2, 1, this, T, result, err ); MSQ_ERRZERO(err);
  wrt_T(2,2) = do_finite_difference( 2, 2, this, T, result, err ); MSQ_ERRZERO(err);
  return true;
}

bool TRel3DMetric::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                       double& result,
                                       MsqMatrix<3,3>& deriv_wrt_T,
                                       MsqMatrix<3,3> hess_wrt_T[6],
                                       MsqError& err )
{
  return do_numerical_hessian( this, T, result, deriv_wrt_T, hess_wrt_T, err );
}

} // namespace Mesquite

