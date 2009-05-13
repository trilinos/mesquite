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


/** \file Target2DShapeSizeOrientBarrierAlt2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeOrientBarrierAlt2.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

msq_std::string Target2DShapeSizeOrientBarrierAlt2::get_name() const
  { return "ShapeSizeOrientBarrier"; }

bool Target2DShapeSizeOrientBarrierAlt2::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError& )
{
  if (invalid_determinant(det(A))) {
    result = 0.0;
    return false;
  }
  MsqMatrix<2,2> T_inv = W * inverse(A);
  T_inv(0,0) -= 1.0;
  T_inv(1,1) -= 1.0;
  result = sqr_Frobenius(T_inv);
  return true;
}


} // namespace Mesquite
