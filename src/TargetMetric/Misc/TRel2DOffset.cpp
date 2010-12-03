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


/** \file TRel2DOffset.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel2DOffset.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

std::string TRel2DOffset::get_name() const
  { return "offset(" + mMetric->get_name() + ')'; }

bool TRel2DOffset::evaluate( const MsqMatrix<2,2>& T, 
                             double& result, 
                             MsqError& err )
{
  bool rval = mMetric->evaluate( T, result, err );
  result += mAlpha;
  return rval;
}

bool TRel2DOffset::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                       double& result,
                                       MsqMatrix<2,2>& deriv_wrt_T,
                                       MsqError& err )
{
  bool rval = mMetric->evaluate_with_grad( T, result, deriv_wrt_T, err );
  result += mAlpha;
  return rval;
}

bool TRel2DOffset::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                       double& result,
                                       MsqMatrix<2,2>& deriv_wrt_T,
                                       MsqMatrix<2,2> second_wrt_T[3],
                                       MsqError& err )
{
  bool rval = mMetric->evaluate_with_hess( T, result, deriv_wrt_T, second_wrt_T, err );
  result += mAlpha;
  return rval;
}


} // namespace MESQUITE_NS
