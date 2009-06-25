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


/** \file Target2DShapeSizeAlt2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeAlt2.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

msq_std::string Target2DShapeSizeAlt2::get_name() const
  { return "ShapeSize"; }

bool Target2DShapeSizeAlt2::evaluate( const MsqMatrix<2,2>& A, 
                                      const MsqMatrix<2,2>& W, 
                                      double& result, 
                                      MsqError&  )
{
  const MsqMatrix<2,2> T = A * inverse(W);
  const double nT = sqr_Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau - 1;
  result = 2*nT - 4*tau + mGamma*tau1*tau1;
  return true;
}


bool Target2DShapeSizeAlt2::evaluate_with_grad( const MsqMatrix<2,2>& A, 
                                                const MsqMatrix<2,2>& W, 
                                                double& result, 
                                                MsqMatrix<2,2>& deriv_wrt_A,
                                                MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double nT = sqr_Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau - 1;
  result = 2*nT - 4*tau + mGamma*tau1*tau1;
 
  deriv_wrt_A = T;
  deriv_wrt_A *= 4;
  deriv_wrt_A += (2*mGamma*tau1 - 4) * transpose_adj(T);
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);

  return true;
}

bool Target2DShapeSizeAlt2::evaluate_with_hess( const MsqMatrix<2,2>& A, 
                                                const MsqMatrix<2,2>& W, 
                                                double& result, 
                                                MsqMatrix<2,2>& deriv_wrt_A,
                                                MsqMatrix<2,2> second[3],
                                                MsqError& err )
{
  const MsqMatrix<2,2> Winv = inverse(W);
  const MsqMatrix<2,2> T = A * Winv;
  const double nT = sqr_Frobenius(T);
  const double tau = det(T);
  const double tau1 = tau - 1;
  result = 2*nT - 4*tau + mGamma*tau1*tau1;
 
  const double f = 2*mGamma*tau1 - 4;
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  deriv_wrt_A = T;
  deriv_wrt_A *= 4;
  deriv_wrt_A += f * adjt;
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);

  set_scaled_outer_product( second, 2*mGamma, adjt );
  pluseq_scaled_I( second, 4 );
  pluseq_scaled_2nd_deriv_of_det( second, f );
  second_deriv_wrt_product_factor( second, Winv );
  
  return true;
}

} // namespace Mesquite
