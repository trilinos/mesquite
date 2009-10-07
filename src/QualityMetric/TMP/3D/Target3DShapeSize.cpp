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
  const MsqMatrix<3,3> T = A * inverse(W);
  const double nT = Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau-1;
  result = nT*nT*nT - 3*MSQ_SQRT_THREE*tau + mGamma*tau1*tau1;
  return true;
}


bool Target3DShapeSize::evaluate_with_grad( const MsqMatrix<3,3>& A, 
                                            const MsqMatrix<3,3>& W, 
                                            double& result, 
                                            MsqMatrix<3,3>& wrt_A,
                                            MsqError&  )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double nT = Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau-1;
  result = nT*nT*nT - 3*MSQ_SQRT_THREE*tau + mGamma*tau1*tau1;
  
  wrt_A = T;
  wrt_A *= 3*nT;
  wrt_A -= (3*MSQ_SQRT_THREE - 2*mGamma*tau1) * transpose_adj(T);
  wrt_A = wrt_A * transpose(Winv);
  
  return true;
}

bool Target3DShapeSize::evaluate_with_hess( const MsqMatrix<3,3>& A,
                                            const MsqMatrix<3,3>& W,
                                            double& result,
                                            MsqMatrix<3,3>& wrt_A,
                                            MsqMatrix<3,3> second[6],
                                            MsqError& err )
{
  const MsqMatrix<3,3> Winv = inverse(W);
  const MsqMatrix<3,3> T = A * Winv;
  const double nT = Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau-1;
  result = nT*nT*nT - 3*MSQ_SQRT_THREE*tau + mGamma*tau1*tau1;
  
  const double f = (3*MSQ_SQRT_THREE - 2*mGamma*tau1);
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  wrt_A = T;
  wrt_A *= 3*nT;
  wrt_A -= f * adjt;
  wrt_A = wrt_A * transpose(Winv);
  
  set_scaled_outer_product( second, 2 * mGamma, adjt );
  pluseq_scaled_2nd_deriv_of_det( second, -f, T );
  pluseq_scaled_I( second, 3*nT );
    // Could perturb T a bit if the norm is zero, but that would just
    // result in the coefficent of the outer product being practically 
    // zero, so just skip the outer product in that case.
    // Anyway nT approaches zero as T does, so the limit of this term
    // as nT approaches zero is zero.
  if (nT > 1e-100)  // NOTE: nT is always positive
    pluseq_scaled_outer_product( second, 3/nT, T );
  
  second_deriv_wrt_product_factor( second, Winv );
  
  return true;
}

} // namespace Mesquite
