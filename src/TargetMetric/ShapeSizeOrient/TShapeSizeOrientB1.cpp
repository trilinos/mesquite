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


/** \file TRel2DShapeSizeOrientBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSizeOrientB1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TShapeSizeOrientB1::get_name() const
  { return "TShapeSizeOrientB1"; }

TShapeSizeOrientB1::~TShapeSizeOrientB1() {}

template <unsigned DIM> static inline
bool eval( const MsqMatrix<DIM,DIM>& T, double& result)
{
  double tau = det(T);
  if (TMetric::invalid_determinant(tau)) {
    result = 0.0;
    return false;
  }
  
  MsqMatrix<DIM,DIM> T_I(T);
  pluseq_scaled_I( T_I, -1 );
  result = sqr_Frobenius( T_I ) / (2*tau);
  return true;
}

template <unsigned DIM> static inline
bool grad( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv )
{
  const double d = det(T);
  if (TMetric::invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  
  deriv = T;
  pluseq_scaled_I( deriv, -1 );
  double inv_d = 1.0/d;
  result = 0.5 * sqr_Frobenius(deriv) * inv_d;
  
  deriv -= result * transpose_adj(T);
  deriv *= inv_d;
  
  return true;
}

template <unsigned DIM> static inline
bool hess( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv, 
           MsqMatrix<DIM,DIM>* second )
{
  const double d = det(T);
  if (TMetric::invalid_determinant(d)) { // barrier
    result = 0.0;
    return false;
  }
  
  deriv = T;
  pluseq_scaled_I( deriv, -1.0 );
  double inv_d = 1.0/d;
  result = 0.5 * sqr_Frobenius(deriv) * inv_d;
  
  MsqMatrix<DIM,DIM> adjt = transpose_adj(T);
  set_scaled_outer_product( second, 2*result*inv_d*inv_d, adjt );
  pluseq_scaled_sum_outer_product( second, -inv_d*inv_d, deriv, adjt );
  pluseq_scaled_2nd_deriv_of_det( second, -result * inv_d, T );
  pluseq_scaled_I( second, inv_d );
  
  deriv -= result * adjt;
  deriv *= inv_d;

  return true;
}

TMP_T_TEMPL_IMPL_COMMON(TShapeSizeOrientB1)

} // namespace Mesquite
