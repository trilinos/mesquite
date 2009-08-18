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


/** \file Target3DShapeBarrierAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DShapeBarrierAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

Target3DShapeBarrierAlt1::~Target3DShapeBarrierAlt1() {}

msq_std::string Target3DShapeBarrierAlt1::get_name() const
  { return "ShapeBarrier1"; }

// \mu_3(T) = \frac{ |T|^2 |adj(T)|^2 } {9 \tau^2} - 1
bool Target3DShapeBarrierAlt1::evaluate( const MsqMatrix<3,3>& A, 
                                         const MsqMatrix<3,3>& W, 
                                         double& result, 
                                         MsqError& )
{
  MsqMatrix<3,3> T = A * inverse(W);
  double f = sqr_Frobenius(T);
  double g = sqr_Frobenius(adj(T));
  double d = det(T);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = (f*g) / (9*d*d) - 1;
  return true;
}


bool Target3DShapeBarrierAlt1::evaluate_with_grad( const MsqMatrix<3,3>& A, 
                                     const MsqMatrix<3,3>& W, 
                                     double& result, 
                                     MsqMatrix<3,3>& wrt_A,
                                     MsqError&  )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  double f = sqr_Frobenius(T);
  double g = sqr_Frobenius(adj(T));
  double d = det(T);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = (f*g) / (9*d*d) - 1;
  
  wrt_A = T;
  wrt_A *= (g + f*f);
  wrt_A -= f * (T * transpose(T) * T);
  wrt_A -= f * g / d * transpose_adj(T);
  wrt_A *= 2 / (9*d*d);
  wrt_A = wrt_A * transpose(Winv);
  
  return true;
}


bool Target3DShapeBarrierAlt1::evaluate_with_hess( const MsqMatrix<3,3>& A,
                                                const MsqMatrix<3,3>& W,
                                                double& result,
                                                MsqMatrix<3,3>& wrt_A,
                                                MsqMatrix<3,3> second[6],
                                                MsqError& err )
{
  MsqMatrix<3,3> Winv = inverse(W);
  MsqMatrix<3,3> T = A * Winv;
  double f = sqr_Frobenius(T);
  double g = sqr_Frobenius(adj(T));
  double d = det(T);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  const double den = 1.0/(9*d*d);
  result = f*g*den- 1;
  
  MsqMatrix<3,3> dg = 2 * (f * T - T * transpose(T) * T);
  MsqMatrix<3,3> df = 2 * T;
  MsqMatrix<3,3> dtau = transpose_adj(T);
  
  wrt_A = g*df + f*dg - 2*f*g/d * transpose_adj(T);
  wrt_A *= den;
  wrt_A = wrt_A * transpose(Winv);
  
  set_scaled_2nd_deriv_norm_sqr_adj( second, den*f, T );
  pluseq_scaled_I( second, 2*den*g );
  pluseq_scaled_sum_outer_product( second, den, dg, df );
  pluseq_scaled_sum_outer_product( second, -2*den*g/d, df, dtau );
  pluseq_scaled_sum_outer_product( second, -2*den*f/d, dg, dtau );
  pluseq_scaled_outer_product( second, 6*den*f*g/(d*d), dtau );
  pluseq_scaled_2nd_deriv_of_det( second, -2*den*f*g/d, T );
  second_deriv_wrt_product_factor( second, Winv );
  
  return true;
}


} // namespace Mesquite
