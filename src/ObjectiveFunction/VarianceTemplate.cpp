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


/** \file VarianceTemplate.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "VarianceTemplate.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

ObjectiveFunction* VarianceTemplate::clone() const
  { return new VarianceTemplate(*this); }
  
void VarianceTemplate::clear()
{
  mCount = 0;
  mSum = mSqrSum = 0;
  saveCount = 0;
  saveSum = saveSqrSum = 0;
}

void VarianceTemplate::accumulate( double sum, 
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

bool VarianceTemplate::evaluate( EvalType type, 
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

bool VarianceTemplate::evaluate_with_gradient( EvalType type, 
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

bool VarianceTemplate::evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        msq_std::vector<Vector3D>& grad_out,
                                        msq_std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient and Hessian data
  grad_out.clear();     // store sum of metric * gradient of metric, and later OF gradient
  grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  tmpGradient.clear();  // store sum of gradients of metrics
  tmpGradient.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  tmpDiag1.clear();     // store sum of Hessians of metrics
  tmpDiag1.resize( pd.num_free_vertices(), SymMatrix3D(0.0) );
  tmpDiag2.clear();     // store sum of metric * Hessian of metric
  tmpDiag2.resize( pd.num_free_vertices(), SymMatrix3D(0.0) );
  hess_diag_out.clear(); // store sum of metric * outer_product(metric gradient), and later OF Hessian
  hess_diag_out.resize( pd.num_free_vertices(), SymMatrix3D(0.0) );
  
    // calculate OF value and gradient for just the patch
  Matrix3D op;
  msq_std::vector<size_t>::const_iterator i;
  double value, sum = 0.0, sqr = 0.0;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_Hessian_diagonal( pd, *i, value, mIndices, mGradient, mHessDiag, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    sum += value;
    sqr += value*value;

    for (size_t j = 0; j < mIndices.size(); ++j) {
      const size_t r = mIndices[j];
      tmpGradient[r] += mGradient[j];
      mGradient[j] *= value;
      grad_out[r] += mGradient[j];
      
      hess_diag_out[r] += outer( mGradient[j] );
      tmpDiag1[r] += mHessDiag[j];
      mHessDiag[j] *= value;
      tmpDiag2[r] += mHessDiag[j];
    }
  }
  
    // update member data
  size_t n;
  accumulate( sum, sqr, qmHandles.size(), type, sum, sqr, n );
  if (n < 2) {
    value_out = 0.0;
    grad_out.clear();
    grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
    hess_diag_out.clear();
    hess_diag_out.resize( pd.num_free_vertices(), SymMatrix3D(0.0) );
    return true;
  }

    // calculate OF value
  value_out = qm->get_negate_flag() * (n*sqr - sum*sum) / (n*(n - 1));
  
    // Finish calculation of gradient and Hessian
  const double dneg = qm->get_negate_flag() * 2.0;
  const double n_inv = 1.0/n;
  const double nless1_inv = 1.0/(n-1);
  const double avg = sum * n_inv;
  const double f = dneg * nless1_inv;
  const double f2 = avg * nless1_inv;
  for (size_t k = 0; k < pd.num_free_vertices(); ++k) {
    tmpGradient[k] *= avg;
    grad_out[k] -= tmpGradient[k];
    grad_out[k] *= f;

    hess_diag_out[k] *= n_inv;
    tmpDiag1[k] *= f2;
    tmpDiag2[k] *= nless1_inv;
    hess_diag_out[k] += tmpDiag1[k];
    hess_diag_out[k] += tmpDiag2[k];
    hess_diag_out[k] *= dneg;
  }
    
  return true;
}


} // namespace Mesquite
