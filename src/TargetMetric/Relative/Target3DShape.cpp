/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009
 Sandia National Laboratories.  Developed at the
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


/** \file Target3DShape.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DShape.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

Target3DShape::~Target3DShape() {}

std::string Target3DShape::get_name() const
  { return "Shape"; }

bool Target3DShape::evaluate( const MsqMatrix<3,3>& A, 
                              const MsqMatrix<3,3>& W, 
                              double& result, 
                              MsqError& )
{
  MsqMatrix<3,3> T = A * inverse(W);
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;
  return true;
}


bool Target3DShape::evaluate_with_grad( const MsqMatrix<3,3>& A, 
                                        const MsqMatrix<3,3>& W, 
                                        double& result, 
                                        MsqMatrix<3,3>& deriv_wrt_A,
                                        MsqError& err )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;

  deriv_wrt_A = T;
  deriv_wrt_A *= f;
  deriv_wrt_A -= MSQ_SQRT_THREE*transpose_adj(T);
  deriv_wrt_A *= 3;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  return true;
}

bool Target3DShape::evaluate_with_hess( const MsqMatrix<3,3>& A, 
                                        const MsqMatrix<3,3>& W, 
                                        double& result, 
                                        MsqMatrix<3,3>& deriv_wrt_A,
                                        MsqMatrix<3,3> second_wrt_A[6],
                                        MsqError& err )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;

  deriv_wrt_A = T;
  deriv_wrt_A *= f;
  deriv_wrt_A -= MSQ_SQRT_THREE*transpose_adj(T);
  deriv_wrt_A *= 3;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  
  set_scaled_2nd_deriv_of_det( second_wrt_A, -3 * MSQ_SQRT_THREE, T );
  pluseq_scaled_outer_product( second_wrt_A, 3.0/f, T );
  pluseq_scaled_I( second_wrt_A, 3.0*f );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );
  return true;
}

} // namespace Mesquite
