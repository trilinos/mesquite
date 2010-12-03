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


/** \file TRel3DSum.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel3DSum.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

std::string TRel3DSum::get_name() const
  { return mu1->get_name() + '+' + mu2->get_name(); }

bool TRel3DSum::evaluate( const MsqMatrix<3,3>& T, 
                            double& result, 
                            MsqError& err )
{
  double val2;
  bool rval = mu1->evaluate( T, result, err );  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate( T, val2, err ); MSQ_ERRZERO(err);
  result += val2;
  return rval && rval2;
}

bool TRel3DSum::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                      double& result,
                                      MsqMatrix<3,3>& deriv_wrt_T,
                                      MsqError& err )
{
  double val2;
  MsqMatrix<3,3> grad2;
  bool rval = mu1->evaluate_with_grad( T, result, deriv_wrt_T, err );  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate_with_grad( T, val2, grad2, err ); MSQ_ERRZERO(err);
  result += val2;
  deriv_wrt_T += grad2;
  return rval && rval2;
}

bool TRel3DSum::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                      double& result,
                                      MsqMatrix<3,3>& deriv_wrt_T,
                                      MsqMatrix<3,3> second_wrt_T[6],
                                      MsqError& err )
{
  double val2;
  MsqMatrix<3,3> grad2, hess2[6];
  bool rval = mu1->evaluate_with_hess( T, result, deriv_wrt_T, second_wrt_T, err );  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate_with_hess( T, val2, grad2, hess2, err ); MSQ_ERRZERO(err);
  result += val2;
  deriv_wrt_T += grad2;
  second_wrt_T[0] += hess2[0];
  second_wrt_T[1] += hess2[1];
  second_wrt_T[2] += hess2[2];
  second_wrt_T[3] += hess2[3];
  second_wrt_T[4] += hess2[4];
  second_wrt_T[5] += hess2[5];
  return rval && rval2;
}


} // namespace MESQUITE_NS
