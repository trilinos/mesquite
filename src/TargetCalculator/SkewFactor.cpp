/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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


/** \file SkewFactor.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "SkewFactor.hpp"
#include "TargetCalculator.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {


bool SkewFactor::get_skew_2D( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<2,2>& M_out,
                              MsqError& err )
{
  MsqMatrix<3,2> W;
  bool result = srcTargets->get_2D_target( pd, element, sample, W, err ); 
  return !MSQ_CHKERR(err) && result && factor_skew( W, M_out );
}

bool SkewFactor::get_skew_3D( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,3>& M_out,
                              MsqError& err )
{
  MsqMatrix<3,3> W;
  bool result = srcTargets->get_3D_target( pd, element, sample, W, err ); 
  return !MSQ_CHKERR(err) && result && factor_skew( W, M_out );
}



} // namespace MESQUITE_NS
