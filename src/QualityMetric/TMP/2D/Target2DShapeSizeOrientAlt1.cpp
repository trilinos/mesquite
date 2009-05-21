/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file Target2DShapeSizeOrientAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeOrientAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

msq_std::string Target2DShapeSizeOrientAlt1::get_name() const
  { return "ShapeSizeOrient"; }

bool Target2DShapeSizeOrientAlt1::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError&  )
{
  result = sqr_Frobenius( A - W );
  return true;
}

bool Target2DShapeSizeOrientAlt1::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                               const MsqMatrix<2,2>& W,
                                               double& result,
                                               MsqMatrix<2,2>& deriv_wrt_A,
                                               MsqError& err )
{
  MsqMatrix<2,2> diff = A - W;
  result = sqr_Frobenius( diff );
  deriv_wrt_A = diff;
  deriv_wrt_A *= 2.0;
  return true;
}

bool Target2DShapeSizeOrientAlt1::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                               const MsqMatrix<2,2>& W,
                                               double& result,
                                               MsqMatrix<2,2>& deriv_wrt_A,
                                               MsqMatrix<2,2> second_wrt_A[3],
                                               MsqError& err )
{
  MsqMatrix<2,2> diff = A - W;
  result = sqr_Frobenius( diff );
  deriv_wrt_A = diff;
  deriv_wrt_A *= 2.0;
  set_scaled_I( second_wrt_A, 2.0 );
  return true;
}

} // namespace Mesquite
