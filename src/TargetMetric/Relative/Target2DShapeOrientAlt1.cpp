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


/** \file Target2DShapeOrientAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeOrientAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string Target2DShapeOrientAlt1::get_name() const
  { return "ShapeOrient"; }

bool Target2DShapeOrientAlt1::evaluate( const MsqMatrix<2,2>& A, 
                                 const MsqMatrix<2,2>& W, 
                                 double& result, 
                                 MsqError&  )
{
  const MsqMatrix<2,2> T = A * inverse(W);
  const double tr = trace(T);
  result = sqr_Frobenius( T ) - 0.5 * tr * fabs(tr);
  return true;
}

bool Target2DShapeOrientAlt1::evaluate_with_grad( const MsqMatrix<2,2>& A, 
                                              const MsqMatrix<2,2>& W, 
                                              double& result, 
                                              MsqMatrix<2,2>& deriv_wrt_A,
                                              MsqError& err )
{
  MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  const double tr = trace(T);
  result = sqr_Frobenius( T ) - 0.5 * tr * fabs(tr);
  deriv_wrt_A = T;
  deriv_wrt_A *= 2;
  deriv_wrt_A(0,0) -= fabs(tr);
  deriv_wrt_A(1,1) -= fabs(tr);
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  return true;
}

bool Target2DShapeOrientAlt1::evaluate_with_Hess( const MsqMatrix<2,2>& A, 
                                              const MsqMatrix<2,2>& W, 
                                              double& result, 
                                              MsqMatrix<2,2>& deriv_wrt_A,
                                              MsqMatrix<2,2> second_wrt_A[3],
                                              MsqError& err )
{
  MsqMatrix<2,2> Winv = inverse(W);
  MsqMatrix<2,2> T = A * Winv;
  const double tr = trace(T);
  result = sqr_Frobenius( T ) - 0.5 * tr * fabs(tr);
  deriv_wrt_A = T;
  deriv_wrt_A *= 2;
  deriv_wrt_A(0,0) -= fabs(tr);
  deriv_wrt_A(1,1) -= fabs(tr);
  deriv_wrt_A = deriv_wrt_A * transpose(Winv);
  set_scaled_I( second_wrt_A, 2.0 );
  pluseq_scaled_outer_product_I_I( second_wrt_A, tr < 0 ? 1 : -1 );
  second_deriv_wrt_product_factor( second_wrt_A, Winv );
  return true;
}



} // namespace Mesquite
