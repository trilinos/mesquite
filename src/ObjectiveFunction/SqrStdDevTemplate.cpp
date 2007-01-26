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


/** \file SqrStdDevTemplate.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "SqrStdDevTemplate.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "PatchData.hpp"

namespace Mesquite {

ObjectiveFunction* SqrStdDevTemplate::clone() const
  { return new SqrStdDevTemplate(*this); }

bool SqrStdDevTemplate::evaluate( EvalType type, 
                               PatchData& pd,
                               double& value_out,
                               bool free,
                               MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, free, err );  MSQ_ERRFALSE(err);
  
    // calculate OF value for just the patch
  msq_std::vector<size_t>::const_iterator i;
  double value, sum = 0.0, sqr = 0.0;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate( pd, *i, value, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    
    sum += value;
    sqr += value*value;
  }
  
    // get overall OF value, update member data, etc.
  size_t n;
  accumulate( sum, sqr, qmHandles.size(), type, sum, sqr, n );
  if (n < 2) {
    value_out = 0.0;
    return true;
  }
  
  value_out = qm->get_negate_flag() * (n*sqr - sum*sum) / (n*(n - 1));
  return true;
}

bool SqrStdDevTemplate::evaluate_with_gradient( EvalType type, 
                                             PatchData& pd,
                                             double& value_out,
                                             msq_std::vector<Vector3D>& grad_out,
                                             MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient
  grad_out.clear();
  grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  tmpGradient.clear();
  tmpGradient.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  
    // calculate OF value and gradient for just the patch
  msq_std::vector<size_t>::const_iterator i;
  double value, sum = 0.0, sqr = 0.0;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_gradient( pd, *i, value, mIndices, mGradient, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    sum += value;
    sqr += value*value;

    for (size_t j = 0; j < mIndices.size(); ++j) {
      tmpGradient[mIndices[j]] += mGradient[j];
      mGradient[j] *= value;
      grad_out[mIndices[j]] += mGradient[j];
    }
  }
  
    // update member data
  size_t n;
  accumulate( sum, sqr, qmHandles.size(), type, sum, sqr, n );
  if (n < 2) {
    value_out = 0.0;
    grad_out.clear();
    grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
    return true;
  }

    // calculate OF value
  value_out = qm->get_negate_flag() * (n*sqr - sum*sum) / (n*(n - 1));
  
    // calculate gradient
  const double avg = sum/n;
  const double f = 2.0 / (n - 1);
  for (size_t k = 0; k < pd.num_free_vertices(); ++k) {
    tmpGradient[k] *= avg;
    grad_out[k] -= tmpGradient[k];
    grad_out[k] *= f;
  }
    
  return true;
}

bool SqrStdDevTemplate::evaluate_with_Hessian( EvalType type, 
                                             PatchData& pd,
                                             double& value_out,
                                             msq_std::vector<Vector3D>& grad_out,
                                             MsqHessian& Hess_out,
                                             MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient
  grad_out.clear();
  grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  tmpGradient.clear();
  tmpGradient.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  Hess_out.zero_out();
  tmpHessian1.initialize(Hess_out);
  tmpHessian1.zero_out();
  tmpHessian2.initialize(Hess_out);
  tmpHessian2.zero_out();
  
    // calculate OF value and gradient for just the patch
  Matrix3D op;
  msq_std::vector<size_t>::const_iterator i;
  double value, sum = 0.0, sqr = 0.0;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_Hessian( pd, *i, value, mIndices, mGradient, mHessian, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    sum += value;
    sqr += value*value;

    size_t h_idx = 0;
    for (size_t j = 0; j < mIndices.size(); ++j) {
      const size_t r = mIndices[j];
      tmpGradient[r] += mGradient[j];
      mGradient[j] *= value;
      grad_out[r] += mGradient[j];
      for (size_t k = j; k < mIndices.size(); ++k) {
        const size_t c = mIndices[j];
        Hess_out.add( r, c, op.outer_product( mGradient[j], mGradient[k] ), err );
        MSQ_ERRZERO(err);
        tmpHessian1.add( r, c, mHessian[h_idx], err );
        MSQ_ERRZERO(err);
        mHessian[h_idx] *= value;
        tmpHessian2.add( r, c, mHessian[h_idx], err );
        MSQ_ERRZERO(err);
        ++h_idx;
      }
    }
  }
  
    // update member data
  size_t n;
  accumulate( sum, sqr, qmHandles.size(), type, sum, sqr, n );
  if (n < 2) {
    value_out = 0.0;
    grad_out.clear();
    grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
    Hess_out.zero_out();
    return true;
  }

    // calculate OF value
  value_out = qm->get_negate_flag() * (n*sqr - sum*sum) / (n*(n - 1));
  
    // calculate gradient
  const double avg = sum/n;
  const double f = qm->get_negate_flag() * 2.0 / (n - 1);
  for (size_t k = 0; k < pd.num_free_vertices(); ++k) {
    tmpGradient[k] *= avg;
    grad_out[k] -= tmpGradient[k];
    grad_out[k] *= f;
  }
  
    // calculate Hessian
  tmpHessian1.scale( -sum );
  tmpHessian2.scale( n );
  tmpHessian1.add( tmpHessian2 );
  tmpHessian1.scale( 1.0 / (n - 1) );
  Hess_out.add( tmpHessian1 );
  Hess_out.scale( qm->get_negate_flag() * 2.0 / n );
    
  return true;
}


} // namespace Mesquite
