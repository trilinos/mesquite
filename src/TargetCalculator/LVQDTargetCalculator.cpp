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


/** \file LVQDTargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LVQDTargetCalculator.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include <algorithm>
#include <assert.h>

namespace MESQUITE_NS {

int LVQDTargetCalculator::add_source( TargetCalculator* source )
{
  if (!source)
    return -1;
  int idx = std::find( uniqueGuides, uniqueGuides+numUniqueGuides, source ) - uniqueGuides;
  if (idx == numUniqueGuides) {
    assert(idx < 4);
    uniqueGuides[idx] = source;
    ++numUniqueGuides;
  }
  return idx;
}

LVQDTargetCalculator::LVQDTargetCalculator( 
                        TargetCalculator* lambda_source,
                        TargetCalculator* V_source,
                        TargetCalculator* Q_source,
                        TargetCalculator* delta_source )
  : numUniqueGuides(0)
{ 
  lambdaIdx = add_source( lambda_source );
  vIdx      = add_source( V_source );
  qIdx      = add_source( Q_source );
  deltaIdx  = add_source( delta_source );
}

LVQDTargetCalculator::~LVQDTargetCalculator() {}

bool LVQDTargetCalculator::get_3D_target( PatchData& pd, 
                                          size_t element,
                                          Sample sample,
                                          MsqMatrix<3,3>& W_out,
                                          MsqError& err )
{
  double lambda[4];
  MsqMatrix<3,3> V[4], Q[4], delta[4], W;
  bool valid;
  
  for (int i = 0; i < numUniqueGuides; ++i) {
    valid = uniqueGuides[i]->get_3D_target( pd, element, sample, W, err );
    if (MSQ_CHKERR(err) || !valid)
      return false;
    valid = factor_3D( W, lambda[i], V[i], Q[i], delta[i], err );
    if (MSQ_CHKERR(err) || !valid)
      return false;
  }
  
  if (vIdx >= 0) {
    W_out = V[vIdx];
    if (lambdaIdx >= 0)
      W_out *= lambda[lambdaIdx];
    if (qIdx >= 0)
      W_out = W_out * Q[qIdx];
    if (deltaIdx >= 0)
      W_out = W_out * delta[deltaIdx];
  }
  else if (qIdx >= 0) {
    W_out = Q[qIdx];
    if (lambdaIdx >= 0)
      W_out *= lambda[lambdaIdx];
    if (deltaIdx >= 0)
      W_out = W_out * delta[deltaIdx];
  }
  else if (deltaIdx >= 0) {
    W_out = delta[deltaIdx];
    if (lambdaIdx >= 0)
      W_out *= lambda[lambdaIdx];
  }
  else if (lambdaIdx >= 0) {
    W_out = MsqMatrix<3,3>(lambda[lambdaIdx]);
  }
  else {
    W_out = MsqMatrix<3,3>(1.0);
  }

  return true;
}

bool LVQDTargetCalculator::get_2D_target( PatchData& pd, 
                                          size_t element,
                                          Sample sample,
                                          MsqMatrix<3,2>& W_out,
                                          MsqError& err )
{
  double lambda[4];
  MsqMatrix<3,2> V[4], W;
  MsqMatrix<2,2> Q[4], delta[4];
  bool valid;
  
  for (int i = 0; i < numUniqueGuides; ++i) {
    valid = uniqueGuides[i]->get_2D_target( pd, element, sample, W, err );
    if (MSQ_CHKERR(err) || !valid)
      return false;
    valid = factor_2D( W, lambda[i], V[i], Q[i], delta[i], err );
    if (MSQ_CHKERR(err) || !valid)
      return false;
  }
  
  if (vIdx >= 0) {
    W_out = V[vIdx];
    if (lambdaIdx >= 0)
      W_out *= lambda[lambdaIdx];
    if (qIdx >= 0)
      W_out = W_out * Q[qIdx];
    if (deltaIdx >= 0)
      W_out = W_out * delta[deltaIdx];
  }
  else if (qIdx >= 0) {
    W_out(0,0) = Q[qIdx](0,0); W_out(0,1) = Q[qIdx](0,1);
    W_out(1,0) = Q[qIdx](1,0); W_out(1,1) = Q[qIdx](1,1);
    W_out(2,0) = 0.0;          W_out(2,1) = 0.0;
    if (lambdaIdx >= 0)
      W_out *= lambda[lambdaIdx];
    if (deltaIdx >= 0)
      W_out = W_out * delta[deltaIdx];
  }
  else if (deltaIdx >= 0) {
    W_out(0,0) = delta[deltaIdx](0,0); W_out(0,1) = delta[deltaIdx](0,1);
    W_out(1,0) = delta[deltaIdx](1,0); W_out(1,1) = delta[deltaIdx](1,1);
    W_out(2,0) = 0.0;                  W_out(2,1) = 0.0;
    if (lambdaIdx >= 0)
      W_out *= lambda[lambdaIdx];
  }
  else if (lambdaIdx >= 0) {
    W_out = MsqMatrix<3,2>(lambda[lambdaIdx]);
  }
  else {
    W_out = MsqMatrix<3,2>(1.0);
  }

  return true;
}


} // namespace Mesquite
