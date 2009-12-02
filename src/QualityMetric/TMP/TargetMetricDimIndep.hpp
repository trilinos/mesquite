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

namespace MESQUITE_NS {

template <typename TargetMetric, unsigned Dim>
static inline double
do_finite_difference( int r, int c, TargetMetric* metric, 
                      MsqMatrix<Dim, Dim> A, 
                      const MsqMatrix<Dim, Dim>& W,
                      double value, MsqError& err )
{
  const double INITAL_STEP = std::max( 1e-6, fabs(1e-9*value) );
  const double init = A(r,c);
  bool valid;
  double diff_value;
  for (double step = INITAL_STEP; step > std::numeric_limits<double>::epsilon(); step *= 0.1) {
    A(r,c) = init + step;
    valid = metric->evaluate( A, W, diff_value, err ); MSQ_ERRZERO(err);
    if (valid)
      return (diff_value - value) / step;
  }
  
    // If we couldn't find a valid step, try stepping in the other
    // direciton
  for (double step = INITAL_STEP; step > std::numeric_limits<double>::epsilon(); step *= 0.1) {
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

template <typename TargetMetric, unsigned Dim>
static inline bool
do_numerical_hessian( TargetMetric* metric, 
                      MsqMatrix<Dim, Dim> A,
                      const MsqMatrix<Dim, Dim>& W,
                      double& value,
                      MsqMatrix<Dim, Dim>& grad, 
                      MsqMatrix<Dim, Dim>* Hess, 
                      MsqError& err )
{
    // zero hessian data
  const int num_block = Dim * (Dim + 1) / 2;
  for (int i = 0; i < num_block; ++i)
    Hess[i].zero();

    // evaluate gradient for input values
  bool valid = metric->evaluate_with_grad( A, W, value, grad, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
    // do finite difference for each term of A
  const double INITAL_STEP = std::max( 1e-6, fabs(1e-9*value) );
  double value2;
  MsqMatrix<Dim,Dim> grad2;
  for (unsigned r = 0; r < Dim; ++r) {  // for each row of A
    for (unsigned c = 0; c < Dim; ++c) {  // for each column of A
      const double in_val = A(r,c);
      double step;
      for (step = INITAL_STEP; step > std::numeric_limits<double>::epsilon(); step *= 0.1) {
        A(r,c) = in_val + step;
        valid = metric->evaluate_with_grad( A, W, value2, grad2, err );  MSQ_ERRZERO(err);
        if (valid)
          break;
      }
      
        // if no valid step size, try step in other direction
      if (!valid) {
        for (step = -INITAL_STEP; step < -std::numeric_limits<double>::epsilon(); step *= 0.1) {
          A(r,c) = in_val + step;
          valid = metric->evaluate_with_grad( A, W, value2, grad2, err );  MSQ_ERRZERO(err);
          if (valid)
            break;
        }
        
          // if still no valid step size, give up.
        if (!valid) {
          MSQ_SETERR(err)("No valid step size for finite difference of 2D target metric.",
                          MsqError::INTERNAL_ERROR);
          return false;
        }
      }
      
        // restore A.
      A(r,c) = in_val;
      
        // add values into result matrix
        // values of grad2, in row-major order, are a single 9-value row of the Hessian
      grad2 -= grad;
      grad2 /= step;
      for (unsigned b = 0; b < r; ++b) {
        const int idx = Dim*b - b*(b+1)/2 + r;
        Hess[idx].add_column( c, transpose( grad2.row(b) ) );
      }
      for (unsigned b = r; b < Dim; ++b) {
        const int idx = Dim*r - r*(r+1)/2 + b;
        Hess[idx].add_row( c, grad2.row(b) );
      }
    } // for (c)
  } // for (r)
  
    // Values in non-diagonal blocks were added twice.
  for (unsigned r = 0, h = 1; r < Dim-1; ++r, ++h)
    for (unsigned c = r + 1; c < Dim; ++c, ++h)
      Hess[h] *= 0.5;
  
  return true;
}

} // namespace Mesquite

#endif
