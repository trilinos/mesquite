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


/** \file TSquared3D.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TSquared3D.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

TSquared3D::~TSquared3D() {}

std::string TSquared3D::get_name() const
  { return "TSquared"; }

bool TSquared3D::evaluate( const MsqMatrix<3,3>& A, 
                           const MsqMatrix<3,3>& W, 
                           double& result, 
                           MsqError& err )
{
  const MsqMatrix<3,3> T = A * inverse(W);
  result = sqr_Frobenius( T );
  return true;
}


bool TSquared3D::evaluate_with_grad( const MsqMatrix<3,3>& A, 
                                     const MsqMatrix<3,3>& W, 
                                     double& result, 
                                     MsqMatrix<3,3>& wrt_A,
                                     MsqError& )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  result = sqr_Frobenius( T );
  wrt_A = 2*T*transpose(Winv);
  return true;
}

bool TSquared3D::evaluate_with_hess( const MsqMatrix<3,3>& A,
                                     const MsqMatrix<3,3>& W,
                                     double& result,
                                     MsqMatrix<3,3>& deriv_wrt_A,
                                     MsqMatrix<3,3> second_wrt_A[6],
                                     MsqError& err )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  MsqMatrix<3,3> V = 2 * transpose(Winv);
  result = sqr_Frobenius( T );
  deriv_wrt_A = T * V;
    // diagonal blocks
  second_wrt_A[0] = second_wrt_A[3] = second_wrt_A[5] = Winv * V;
    // non-diagonal blocks are zero
  second_wrt_A[1] = second_wrt_A[2] = second_wrt_A[4] = MsqMatrix<3,3>(0.0);
  return true;
}


} // namespace Mesquite
