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


/** \file Target3DShapeSize.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target3DShapeSize.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

Target3DShapeSize::~Target3DShapeSize() {}

std::string Target3DShapeSize::get_name() const
  { return "ShapeSize"; }

bool Target3DShapeSize::evaluate( const MsqMatrix<3,3>& A, 
                                  const MsqMatrix<3,3>& W, 
                                  double& result, 
                                  MsqError& )
{
  MsqMatrix<3,3> T = A * inverse(W);
  double f = Frobenius(T);
  double d = det(T);
  double d1 = d-1;
  double den = 3 * MSQ_SQRT_THREE * d;
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = (f*f*f)/den + mGamma * d1 * d1;
  return true;
}


/** \f$ \frac{\partial}{\partial T} 
 *     = \frac{\left| T \right|}{\sqrt{3} det(T)}T
 *     + \left[2 \gamma \left(det(T) - 1\right) 
 *     - \frac{\left|T\right|^3}{3 \sqrt{3} det^2(T)}\right] adj(T^t) \f$
 */
bool Target3DShapeSize::evaluate_with_grad( const MsqMatrix<3,3>& A, 
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
  double d1 = d-1;
  double den = 3 * MSQ_SQRT_THREE * d;
  double norm_cube = norm*norm*norm;
  result = norm_cube/den + mGamma*d1*d1;
  
  double f1 = norm / (d*MSQ_SQRT_THREE);
  double f2 = 2*mGamma*d1 - norm_cube / (3*MSQ_SQRT_THREE * d*d);
  wrt_A = (f1*T + f2*transpose_adj(T)) * transpose(Winv);
  return true;
}

/**
 *\f$ \frac{\left|T\right|}{\sqrt{3} \tau} I_9 + 
       \left( 2 \gamma (\tau - 1) - \frac{\left|T\right|^3}{3\sqrt{3} \tau^3} \right) 
          \frac{\partial^2 \tau}{\partial T^2} +
       2 \left( \gamma + \frac{\left|T\right|^3}{3 \sqrt{3} \tau^3} \right)
           \left( \frac{\partial \tau}{\partial T} \otimes \frac{\partial \tau}{\partial T} \right) -
         \frac{\left|T\right|}{\sqrt{3} \tau^2}\left( T \otimes \frac{\partial \tau}{\partial T} + 
                \frac{\partial \tau}{\partial T} \otimes T \right) \f$
 */
bool Target3DShapeSize::evaluate_with_hess( const MsqMatrix<3,3>& A,
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
  
  double inv_det = 1/d;
  double den = 1.0 / (MSQ_SQRT_THREE * d);
  double third = 1.0 / 3.0;
  double norm = Frobenius(T);
  double d1 = d-1;
  double norm_cube = norm*norm*norm;
  result = norm_cube*third*den + mGamma*d1*d1;
  
  double f1 = norm * den;
  double f2 = 2 * mGamma * d1 - norm_cube * third * den * inv_det;
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  deriv_wrt_A = (f1*T + f2*adjt) * transpose(Winv);
  
  set_scaled_outer_product( second_wrt_A, den / norm, T );
  pluseq_scaled_I( second_wrt_A, norm * den );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_A, 2 * mGamma * d1 - norm_cube * third * den * inv_det, T );
  pluseq_scaled_outer_product( second_wrt_A, 2*(norm_cube * third * den * inv_det * inv_det + mGamma), adjt );
  pluseq_scaled_sum_outer_product( second_wrt_A, -norm * den * inv_det, T, adjt );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );

  return true;
}

} // namespace Mesquite
