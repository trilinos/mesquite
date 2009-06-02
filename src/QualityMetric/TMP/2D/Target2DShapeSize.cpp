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


/** \file Target2DShapeSize.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSize.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

msq_std::string Target2DShapeSize::get_name() const
  { return "ShapeSize"; }

/** \f$ \mu(T) = |T|^2 - 2 \psi(T) + 2 \f$
 *  \f$ \psi(T) = \sqrt{|T|^2 + 2 \tau \f$
 *  \f$ \tau = det(T) \f$
 */
bool Target2DShapeSize::evaluate( const MsqMatrix<2,2>& A, 
                                  const MsqMatrix<2,2>& W, 
                                  double& result, 
                                  MsqError&  )
{
  const MsqMatrix<2,2> T = A * inverse(W);
  const double frob_sqr = sqr_Frobenius(T);
  result = frob_sqr - 2.0 * sqrt( frob_sqr + 2.0*det(T) ) + 2.0;
  return true;
}


bool Target2DShapeSize::evaluate_with_grad( const MsqMatrix<2,2>& A, 
                                            const MsqMatrix<2,2>& W, 
                                            double& result, 
                                            MsqMatrix<2,2>& deriv_wrt_A,
                                            MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double frob_sqr = sqr_Frobenius(T);
  const double psi = sqrt( frob_sqr + 2.0*det(T) );
  result = frob_sqr - 2.0 * psi + 2.0;

  if (fabs(psi) > DBL_EPSILON) {
    deriv_wrt_A = T;
    deriv_wrt_A *= (1.0 - 1.0/psi);
    deriv_wrt_A -= 1.0/psi * transpose_adj(T);
    deriv_wrt_A *= 2;
    deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  }
  else {
    deriv_wrt_A = MsqMatrix<2,2>(0.0);
  }
  return true;
}

bool Target2DShapeSize::evaluate_with_hess( const MsqMatrix<2,2>& A, 
                                            const MsqMatrix<2,2>& W, 
                                            double& result, 
                                            MsqMatrix<2,2>& deriv_wrt_A,
                                            MsqMatrix<2,2> second[3],
                                            MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double frob_sqr = sqr_Frobenius(T);
  const double psi = sqrt( frob_sqr + 2.0*det(T) );
  result = frob_sqr - 2.0 * psi + 2.0;

  if (fabs(psi) > DBL_EPSILON) {
    deriv_wrt_A = T;
    deriv_wrt_A *= (1.0 - 1.0/psi);
    deriv_wrt_A -= 1.0/psi * transpose_adj(T);
    deriv_wrt_A *= 2;
    deriv_wrt_A = deriv_wrt_A * transpose(Winv);

    const double s = T(0,1) - T(1,0);
    const double tr = trace(T);
    const double f = -2.0/(psi*psi*psi);
    second[0](0,0) = second[1](0,1) = second[2](1,1) =  f*s*s;
    second[0](0,1) = second[0](1,0) = second[1](1,1) = -f*s*tr;
    second[1](0,0) = second[2](0,1) = second[2](1,0) =  f*s*tr;
    second[0](1,1) = second[2](0,0) = -(second[1](1,0) = -f*tr*tr);
    pluseq_scaled_I( second, 2 );
    second_deriv_wrt_product_factor( second, Winv );
  }
  else {
    deriv_wrt_A = second[0] = second[1] = second[2] = MsqMatrix<2,2>(0.0);
  }
  
  return true;
}

} // namespace Mesquite
