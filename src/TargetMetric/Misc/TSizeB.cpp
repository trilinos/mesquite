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


/** \file TSizeBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TSizeB.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TSizeB::get_name() const
  { return "SizeB"; }

template <unsigned DIM> static inline
bool eval( const MsqMatrix<DIM,DIM>& T, double& result )
{
  double d = det(T);
  if (TMetric::invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = d + 1.0/d - 2.0;
  return true;  
}

template <unsigned DIM> static inline
bool grad( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv_wrt_T )
{
  double d = det(T);
  if (TMetric::invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = d + 1.0/d - 2.0;
  deriv_wrt_T = (1 - 1/(d*d)) * transpose_adj(T);
  return true;  
}

template <unsigned DIM> static inline
bool hess( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv_wrt_T, 
           MsqMatrix<DIM,DIM>* second_wrt_T )
{
  double d = det(T);
  if (TMetric::invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = d + 1.0/d - 2.0;
  MsqMatrix<DIM,DIM> adjt = transpose_adj(T);
  const double f = 1 - 1/(d*d);
  deriv_wrt_T = f * adjt;
  
  set_scaled_outer_product( second_wrt_T, 2/(d*d*d), adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, f, T );
  return true;  
}

TMP_T_TEMPL_IMPL_COMMON(TSizeB)

} // namespace Mesquite
