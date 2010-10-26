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


/** \file Target2DShapeSizeOrient.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeOrient.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

std::string Target2DShapeSizeOrient::get_name() const
  { return "ShapeSizeOrient"; }

bool Target2DShapeSizeOrient::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError&  )
{
  MsqMatrix<2,2> T = A * inverse(W);
  T(0,0) -= 1.0;
  T(1,1) -= 1.0;
  result = sqr_Frobenius( T );
  return true;
}

bool Target2DShapeSizeOrient::evaluate_with_grad( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqMatrix<2,2>& wrt_A,
                                 MsqError&  )
{
  MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  T(0,0) -= 1.0;
  T(1,1) -= 1.0;
  result = sqr_Frobenius( T );
  wrt_A = 2 * T * transpose(Winv);
  return true;
}

bool Target2DShapeSizeOrient::evaluate_with_hess( const MsqMatrix<2,2>& A,
                                     const MsqMatrix<2,2>& W,
                                     double& result,
                                     MsqMatrix<2,2>& deriv_wrt_A,
                                     MsqMatrix<2,2> second_wrt_A[3],
                                     MsqError& err )
{
  MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  MsqMatrix<2,2> V = 2 * transpose(Winv);
  T(0,0) -= 1.0;
  T(1,1) -= 1.0;
  result = sqr_Frobenius( T );
  deriv_wrt_A = T * V;
    // diagonal blocks
  second_wrt_A[0] = second_wrt_A[2] = Winv * V;
    // non-diagonal blocks are zero
  second_wrt_A[1].zero();
  return true;
}



} // namespace Mesquite
