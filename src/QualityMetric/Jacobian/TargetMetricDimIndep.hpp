/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
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
 
    (2008) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TargetMetricDimIndep.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef TARGET_METRIC_DIM_INTEP_HPP
#define TARGET_METRIC_DIM_INTEP_HPP

#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include <limits>

namespace Mesquite {

template <typename TargetMetric>
static inline double
do_finite_difference( int r, int c, TargetMetric* metric, 
         MsqMatrix<TargetMetric::MATRIX_DIM, TargetMetric::MATRIX_DIM> A, 
         const MsqMatrix<TargetMetric::MATRIX_DIM, TargetMetric::MATRIX_DIM>& W,
         double value, MsqError& err )
{
  const double init = A(r,c);
  bool valid;
  double diff_value;
  for (double step = 1e-6; step < std::numeric_limits<double>::epsilon(); step *= 0.1) {
    A(r,c) = init + step;
    valid = metric->evaluate( A, W, diff_value, err ); MSQ_ERRZERO(err);
    if (valid)
      return (diff_value - value) / step;
  }
  
    // If we couldn't find a valid step, try stepping in the other
    // direciton
  for (double step = 1e-6; step > std::numeric_limits<double>::epsilon(); step *= 0.1) {
    A(r,c) = init - step;
    valid = metric->evaluate( A, W, diff_value, err ); MSQ_ERRZERO(err);
    if (valid)
      return (value - diff_value) / step;
  }
  
    // If that didn't work either, then give up.
  MSQ_SETERR(err)("No valid step size for finite difference of 2D target metric.",
                  MsqError::INTERNAL_ERROR);
  return 0.0;
}

} // namespace Mesquite

#endif
