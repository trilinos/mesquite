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


/** \file TRel3DShapeBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel3DShapeBarrier.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

TRel3DShapeBarrier::~TRel3DShapeBarrier() {}

std::string TRel3DShapeBarrier::get_name() const
  { return "ShapeBarrier"; }

bool TRel3DShapeBarrier::evaluate( const MsqMatrix<3,3>& T, 
                                   double& result, 
                                   MsqError& )
{
  double f = Frobenius(T);
  double d = det(T);
  double den = 3 * MSQ_SQRT_THREE * d;
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = (f*f*f)/den - 1.0;
  return true;
}


bool TRel3DShapeBarrier::evaluate_with_grad( const MsqMatrix<3,3>& T, 
                                             double& result, 
                                             MsqMatrix<3,3>& wrt_T,
                                             MsqError&  )
{
  double d = det(T);
  if (d < 1e-12)
    return false;
    
  double norm = Frobenius(T);
  double den = 1.0/(3 * MSQ_SQRT_THREE * d);
  double norm_cube = norm*norm*norm;
  result = norm_cube * den - 1.0;
  wrt_T = T;
  wrt_T *= 3 * norm * den;
  wrt_T -= norm_cube * den/d * transpose_adj(T);
   return true;
}

bool TRel3DShapeBarrier::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                             double& result,
                                             MsqMatrix<3,3>& deriv_wrt_T,
                                             MsqMatrix<3,3> second_wrt_T[6],
                                             MsqError& err )
{
  double d = det(T);
  if (d < 1e-12)
    return false;
  
  double id = 1.0/d;
  double norm = Frobenius(T);
  double den = 1.0/(3 * MSQ_SQRT_THREE * d);
  double norm_cube = norm*norm*norm;
  result = norm_cube * den - 1.0;
  MsqMatrix<3,3> adjt = transpose_adj(T);
  deriv_wrt_T = T;
  deriv_wrt_T *= 3 * norm * den;
  deriv_wrt_T -= norm_cube * den * id * transpose_adj(T);
 
  set_scaled_outer_product( second_wrt_T, 3 * den / norm, T );
  pluseq_scaled_I( second_wrt_T, 3 * norm * den );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -den * norm_cube * id, T );
  pluseq_scaled_outer_product( second_wrt_T, 2 * den * norm_cube * id * id , adjt );
  pluseq_scaled_sum_outer_product( second_wrt_T, -3 * norm * den * id, T, adjt );

  return true;
}

} // namespace Mesquite
