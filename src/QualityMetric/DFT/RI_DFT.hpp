/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */

/*! \file RI_DFT.hpp

Header file for the Mesquite::RI_DFT class

  \author Thomas Leurent
  \date   2004-09-29
 */


#ifndef RI_DFT_hpp
#define RI_DFT_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MsqMessage.hpp"
#include "DistanceFromTarget.hpp"

namespace Mesquite
{
  
  /*! \class RI_DFT
    \brief Class containing the target corner matrices for the context based smoothing. 
  */
  class RI_DFT : public DistanceFromTarget
  {
  public:

    RI_DFT()
    {
      MsqError err;
      set_averaging_method(LINEAR, err); MSQ_CHKERR(err);
      set_metric_type(ELEMENT_BASED);
      set_gradient_type(NUMERICAL_GRADIENT);
      set_hessian_type(NUMERICAL_HESSIAN);
    }
    
    //! virtual destructor ensures use of polymorphism during destruction
    virtual ~RI_DFT()
       {};

    virtual bool evaluate_element(PatchData& pd,
                            MsqMeshEntity* element,
                            double& value, MsqError &err);
    
  protected:
 
  private:
    
  };

  
} //namespace


#endif // RI_DFT_hpp
