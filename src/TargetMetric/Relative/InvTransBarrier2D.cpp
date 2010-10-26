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


/** \file InvTransBarrier2D.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "InvTransBarrier2D.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

std::string InvTransBarrier2D::get_name() const
  { return "InvTransBarrier"; }

bool InvTransBarrier2D::evaluate( const MsqMatrix<2,2>& A, 
                                  const MsqMatrix<2,2>& W, 
                                  double& result, MsqError& err )
{
  double da = det(A);
  double dw = det(W);
  if (invalid_determinant(da) || invalid_determinant(dw))
    return false;
  MsqMatrix<2,2> Ap = transpose_adj(A);
  Ap *= 1.0/da;
  MsqMatrix<2,2> Wp = transpose_adj(W);
  Wp *= 1.0/dw;
  bool rval = metricPtr->evaluate( Ap, Wp, result, err );
  return !MSQ_CHKERR(err) && rval;
}


} // namespace Mesquite
