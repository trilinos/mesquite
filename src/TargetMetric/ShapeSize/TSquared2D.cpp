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


/** \file TSquared2D.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TSquared2D.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

std::string TSquared2D::get_name() const
  { return "TSquared"; }

bool TSquared2D::evaluate( const MsqMatrix<2,2>& T, 
                          double& result, MsqError& )
{
  result = sqr_Frobenius( T );
  return true;
}

bool TSquared2D::evaluate_with_grad( const MsqMatrix<2,2>& T, 
                                     double& result, 
                                     MsqMatrix<2,2>& wrt_T,
                                     MsqError& )
{
  result = sqr_Frobenius( T );
  wrt_T = 2*T;
  return true;
}

bool TSquared2D::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                     double& result,
                                     MsqMatrix<2,2>& deriv_wrt_T,
                                     MsqMatrix<2,2> second_wrt_T[3],
                                     MsqError& err )
{
  result = sqr_Frobenius( T );
  deriv_wrt_T = 2 * T;
    // diagonal blocks
  second_wrt_T[0] = second_wrt_T[2] = MsqMatrix<2,2>(2.0);
    // non-diagonal blocks are zero
  second_wrt_T[1].zero();
  return true;
}


} // namespace Mesquite
