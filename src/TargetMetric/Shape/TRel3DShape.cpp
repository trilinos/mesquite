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


/** \file TRel3DShape.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel3DShape.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

TRel3DShape::~TRel3DShape() {}

std::string TRel3DShape::get_name() const
  { return "Shape"; }

bool TRel3DShape::evaluate( const MsqMatrix<3,3>& T, 
                            double& result, 
                            MsqError& )
{
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;
  return true;
}


bool TRel3DShape::evaluate_with_grad( const MsqMatrix<3,3>& T, 
                                      double& result, 
                                      MsqMatrix<3,3>& deriv_wrt_T,
                                      MsqError& err )
{
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;

  deriv_wrt_T = T;
  deriv_wrt_T *= f;
  deriv_wrt_T -= MSQ_SQRT_THREE*transpose_adj(T);
  deriv_wrt_T *= 3;
  return true;
}

bool TRel3DShape::evaluate_with_hess( const MsqMatrix<3,3>& T, 
                                      double& result, 
                                      MsqMatrix<3,3>& deriv_wrt_T,
                                      MsqMatrix<3,3> second_wrt_T[6],
                                      MsqError& err )
{
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;

  deriv_wrt_T = T;
  deriv_wrt_T *= f;
  deriv_wrt_T -= MSQ_SQRT_THREE*transpose_adj(T);
  deriv_wrt_T *= 3;
  
  set_scaled_2nd_deriv_of_det( second_wrt_T, -3 * MSQ_SQRT_THREE, T );
  pluseq_scaled_outer_product( second_wrt_T, 3.0/f, T );
  pluseq_scaled_I( second_wrt_T, 3.0*f );
  return true;
}

} // namespace Mesquite
