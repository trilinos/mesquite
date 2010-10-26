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


/** \file Target3DSum.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DSum.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

std::string Target3DSum::get_name() const
  { return mu1->get_name() + '+' + mu2->get_name(); }

bool Target3DSum::evaluate( const MsqMatrix<3,3>& A, 
                            const MsqMatrix<3,3>& W, 
                            double& result, 
                            MsqError& err )
{
  double val2;
  bool rval = mu1->evaluate( A, W, result, err );  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate( A, W, val2, err ); MSQ_ERRZERO(err);
  result += val2;
  return rval && rval2;
}

bool Target3DSum::evaluate_with_grad( const MsqMatrix<3,3>& A,
                                      const MsqMatrix<3,3>& W,
                                      double& result,
                                      MsqMatrix<3,3>& deriv_wrt_A,
                                      MsqError& err )
{
  double val2;
  MsqMatrix<3,3> grad2;
  bool rval = mu1->evaluate_with_grad( A, W, result, deriv_wrt_A, err );  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate_with_grad( A, W, val2, grad2, err ); MSQ_ERRZERO(err);
  result += val2;
  deriv_wrt_A += grad2;
  return rval && rval2;
}

bool Target3DSum::evaluate_with_hess( const MsqMatrix<3,3>& A,
                                      const MsqMatrix<3,3>& W,
                                      double& result,
                                      MsqMatrix<3,3>& deriv_wrt_A,
                                      MsqMatrix<3,3> second_wrt_A[6],
                                      MsqError& err )
{
  double val2;
  MsqMatrix<3,3> grad2, hess2[6];
  bool rval = mu1->evaluate_with_hess( A, W, result, deriv_wrt_A, second_wrt_A, err );  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate_with_hess( A, W, val2, grad2, hess2, err ); MSQ_ERRZERO(err);
  result += val2;
  deriv_wrt_A += grad2;
  second_wrt_A[0] += hess2[0];
  second_wrt_A[1] += hess2[1];
  second_wrt_A[2] += hess2[2];
  second_wrt_A[3] += hess2[3];
  second_wrt_A[4] += hess2[4];
  second_wrt_A[5] += hess2[5];
  return rval && rval2;
}


} // namespace MESQUITE_NS
