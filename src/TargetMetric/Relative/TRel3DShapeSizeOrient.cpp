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


/** \file TRel3DShapeSizeOrient.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel3DShapeSizeOrient.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

TRel3DShapeSizeOrient::~TRel3DShapeSizeOrient() {}

std::string TRel3DShapeSizeOrient::get_name() const
  { return "ShapeSizeOrient"; }

bool TRel3DShapeSizeOrient::evaluate( const MsqMatrix<3,3>& T, 
                                        double& result, 
                                        MsqError& err )
{
  MsqMatrix<3,3> T_I(T);
  T_I(0,0) -= 1.0;
  T_I(1,1) -= 1.0;
  T_I(2,2) -= 1.0;
  result = sqr_Frobenius( T_I );
  return true;
}


bool TRel3DShapeSizeOrient::evaluate_with_grad( const MsqMatrix<3,3>& T, 
                                 double& result, 
                                 MsqMatrix<3,3>& wrt_T,
                                 MsqError&  )
{
  MsqMatrix<3,3> T_I(T);
  T_I(0,0) -= 1.0;
  T_I(1,1) -= 1.0;
  T_I(2,2) -= 1.0;
  result = sqr_Frobenius( T_I );
  wrt_T = 2 * T_I;
  return true;
}

bool TRel3DShapeSizeOrient::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                     double& result,
                                     MsqMatrix<3,3>& deriv_wrt_T,
                                     MsqMatrix<3,3> second_wrt_T[6],
                                     MsqError& err )
{
  MsqMatrix<3,3> T_I(T);
  T_I(0,0) -= 1.0;
  T_I(1,1) -= 1.0;
  T_I(2,2) -= 1.0;
  result = sqr_Frobenius( T_I );
  deriv_wrt_T = 2 * T_I;
    // diagonal blocks
  second_wrt_T[0] = second_wrt_T[3] = second_wrt_T[5] = MsqMatrix<3,3>(2.0);
    // non-diagonal blocks are zero
  second_wrt_T[1] = second_wrt_T[2] = second_wrt_T[4] = MsqMatrix<3,3>(0.0);
  return true;
}

} // namespace Mesquite
