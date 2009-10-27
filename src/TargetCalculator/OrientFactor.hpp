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


/** \file OrientFactor.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_ORIENT_FACTOR_HPP
#define MSQ_ORIENT_FACTOR_HPP

#include "Mesquite.hpp"
#include "TargetOrientation.hpp"

namespace MESQUITE_NS {

class TargetCalculator;

/**\brief Calculate orientation component of target matrix from some
 *        other target matrix.
 */
class OrientFactor : public TargetOrientation {
  private:
    TargetCalculator* srcTargets;
  
  public:
    OrientFactor( TargetCalculator* src ) : srcTargets(src) {}
  
    virtual 
    bool get_orient_2D( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,2>& M_out,
                        MsqError& err );
  
    virtual 
    bool get_orient_3D( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,3>& M_out,
                        MsqError& err );
};


} // namespace MESQUITE_NS

#endif
