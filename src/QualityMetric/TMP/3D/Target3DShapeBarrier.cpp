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


/** \file Target3DShapeBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DShapeBarrier.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

Target3DShapeBarrier::~Target3DShapeBarrier() {}

msq_std::string Target3DShapeBarrier::get_name() const
  { return "ShapeBarrier"; }

bool Target3DShapeBarrier::evaluate( const MsqMatrix<3,3>& A, 
                                     const MsqMatrix<3,3>& W, 
                                     double& result, 
                                     MsqError& )
{
  MsqMatrix<3,3> T = A * inverse(W);
  double f = Frobenius(T);
  double d = det(T);
  double den = 3 * MSQ_SQRT_THREE * d;
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = (f*f*f)/den;
  return true;
}


bool Target3DShapeBarrier::evaluate_with_grad( const MsqMatrix<3,3>& A, 
                                 const MsqMatrix<3,3>& W, 
                                 double& result, 
                                 MsqMatrix<3,3>& wrt_A,
                                 MsqError&  )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  double d = det(T);
  if (d < 1e-12)
    return false;
    
  double norm = Frobenius(T);
  double den = 1.0/(3 * MSQ_SQRT_THREE * d);
  double norm_cube = norm*norm*norm;
  result = norm_cube * den;
  wrt_A = T;
  wrt_A *= 3 * norm * den;
  wrt_A -= norm_cube * den/d * transpose_adj(T);
  wrt_A = wrt_A * transpose(Winv);
  return true;
}

bool Target3DShapeBarrier::evaluate_with_hess( const MsqMatrix<3,3>& A,
                                            const MsqMatrix<3,3>& W,
                                            double& result,
                                            MsqMatrix<3,3>& deriv_wrt_A,
                                            MsqMatrix<3,3> second_wrt_A[6],
                                            MsqError& err )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  double d = det(T);
  if (d < 1e-12)
    return false;
  
  double id = 1.0/d;
  double norm = Frobenius(T);
  double den = 1.0/(3 * MSQ_SQRT_THREE * d);
  double norm_cube = norm*norm*norm;
  result = norm_cube * den;
  MsqMatrix<3,3> adjt = transpose_adj(T);
  deriv_wrt_A = T;
  deriv_wrt_A *= 3 * norm * den;
  deriv_wrt_A -= norm_cube * den * id * transpose_adj(T);
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
 
  set_scaled_outer_product( second_wrt_A, 3 * den / norm, T );
  pluseq_scaled_I( second_wrt_A, 3 * norm * den );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, -den * norm_cube * id, T );
  pluseq_scaled_outer_product( second_wrt_A, 2 * den * norm_cube * id * id , adjt );
  pluseq_scaled_sum_outer_product( second_wrt_A, -3 * norm * den * id, T, adjt );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );

  return true;
}

} // namespace Mesquite
