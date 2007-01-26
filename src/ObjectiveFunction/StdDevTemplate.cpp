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


/** \file StdDevTemplate.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "StdDevTemplate.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "PatchData.hpp"

namespace Mesquite {

ObjectiveFunction* StdDevTemplate::clone() const
  { return new StdDevTemplate(*this); }
  
void StdDevTemplate::clear()
{
  mCount = 0;
  mSum = mSqrSum = 0;
  saveCount = 0;
  saveSum = saveSqrSum = 0;
}

void StdDevTemplate::accumulate( double sum, 
                           double sqr_sum,
                           size_t count, 
                           EvalType type,
                           double& result_sum,
                           double& result_sqr,
                           size_t& global_count )
{
  switch (type) 
  {
    case CALCULATE:
      result_sum = sum;
      result_sqr = sqr_sum;
      global_count = count;
      break;
    
    case ACCUMULATE:
      result_sum = mSum += sum;
      result_sqr = mSqrSum += sqr_sum;
      global_count = mCount += count;
      break;
    
    case SAVE:
      saveSum = sum;
      saveSqrSum = sqr_sum;
      saveCount = count;
      result_sum = mSum;
      result_sqr = mSqrSum;
      global_count = mCount;
      break;
    
    case UPDATE:
      mSum -= saveSum;
      mSqrSum -= saveSqrSum;
      mCount -= saveCount;
      result_sum = mSum += saveSum = sum;
      result_sqr = mSqrSum += saveSqrSum = sqr_sum;
      global_count = mCount += saveCount = count;
      break;
    
    case TEMPORARY:
      result_sum = mSum - saveSum + sum;
      result_sqr = mSqrSum - saveSqrSum + sqr_sum;
      global_count = mCount + count - saveCount;
      break;
  }
}

bool StdDevTemplate::evaluate( EvalType type, 
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
  
  value_out = qm->get_negate_flag() * sqrt( (n*sqr - sum*sum)/(n*(n - 1)) );
  return true;
}

bool StdDevTemplate::evaluate_with_gradient( EvalType type, 
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
  const double inv_sqrt_n = 1.0 / sqrt( (double)(n * (n - 1)) );
  const double sqrt_sum = qm->get_negate_flag() * sqrt( n * sqr - sum * sum );
  value_out = sqrt_sum * inv_sqrt_n;
  
    // calculate gradient
  const double f = inv_sqrt_n / sqrt_sum;
  for (size_t k = 0; k < pd.num_free_vertices(); ++k) {
    grad_out[k] *= n;
    tmpGradient[k] *= sum;
    grad_out[k] -= tmpGradient[k];
    grad_out[k] *= f;
  }
    
  return true;
}


} // namespace Mesquite
