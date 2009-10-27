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
#include "TargetSize.hpp"
#include "TargetOrientation.hpp"
#include "TargetSkew.hpp"
#include "TargetAspect.hpp"

namespace MESQUITE_NS {

LVQDTargetCalculator::LVQDTargetCalculator( 
                        TargetSize*        lambda_source,
                        TargetOrientation* V_source,
                        TargetSkew*        Q_source,
                        TargetAspect*      D_source )
  : lambdaGuide( lambda_source ),
    vGuide( V_source ),
    qGuide( Q_source ),
    dGuide( D_source )
{ }

LVQDTargetCalculator::~LVQDTargetCalculator() {}

bool LVQDTargetCalculator::get_3D_target( PatchData& pd, 
                                          size_t element,
                                          Sample sample,
                                          MsqMatrix<3,3>& W_out,
                                          MsqError& err )
{
  double size;
  bool valid;
  if (lambdaGuide) {
    valid = lambdaGuide->get_size( pd, element, sample, size, err ); 
    if (MSQ_CHKERR(err) || !valid)
      return false;
  }
  else {
    size = 1.0;
  }
  W_out = MsqMatrix<3,3>(size);

  MsqMatrix<3,3> M;
  if (vGuide) {
    valid = vGuide->get_orient_3D( pd, element, sample, M, err ); 
    if (MSQ_CHKERR(err) || !valid)
      return false;
    W_out = W_out * M;
  }
  if (qGuide) {
    valid = qGuide->get_skew_3D( pd, element, sample, M, err ); 
    if (MSQ_CHKERR(err) || !valid)
      return false;
    W_out = W_out * M;
  }
  if (dGuide) {
    valid = dGuide->get_aspect_3D( pd, element, sample, M, err ); 
    if (MSQ_CHKERR(err) || !valid)
      return false;
    W_out = W_out * M;
  }

  return true;
}

bool LVQDTargetCalculator::get_2D_target( PatchData& pd, 
                                          size_t element,
                                          Sample sample,
                                          MsqMatrix<3,2>& W_out,
                                          MsqError& err )
{
  double size;
  bool valid;
  if (lambdaGuide) {
    valid = lambdaGuide->get_size( pd, element, sample, size, err ); 
    if (MSQ_CHKERR(err) || !valid)
      return false;
  }
  else {
    size = 1.0;
  }
  MsqMatrix<2,2> M2(size);
  
  if (qGuide) {
    MsqMatrix<2,2> M;
    valid = qGuide->get_skew_2D( pd, element, sample, M, err );
    if (MSQ_CHKERR(err) || !valid)
      return false;
    M2 = M2 * M;
  }
  if (dGuide) {
    MsqMatrix<2,2> M;
    valid = dGuide->get_aspect_2D( pd, element, sample, M, err );
    if (MSQ_CHKERR(err) || !valid)
      return false;
    M2 = M2 * M;
  }
  if (vGuide) {
    valid = vGuide->get_orient_2D( pd, element, sample, W_out, err );
    if (MSQ_CHKERR(err) || !valid)
      return false;
    W_out = W_out * M2;
  }
  else {
    W_out.set_row( 0, M2.row(0) );
    W_out.set_row( 1, M2.row(1) );
    W_out(2,0) = W_out(2,1) = 0.0;
  }

  return true;
}


} // namespace Mesquite
