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

/*! \file DefaultTargetCalculator.hpp

Header file for the Mesquite::DefaultTargetCalculator class

  \author Thomas Leurent
  \date   2004-09-31
 */


#ifndef DefaultTargetCalculator_hpp
#define DefaultTargetCalculator_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MsqMessage.hpp"
#include "TargetCalculator.hpp"

namespace Mesquite
{
  
  /*! \class DefaultTargetCalculator
    \brief This class just compute the default "isotropic" target matrices.
           In fact it delegates this computation to the base class
           function TargetCalculator::compute_default_target_matrices . 
  */
  class DefaultTargetCalculator : public TargetCalculator
  {
  public:
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~DefaultTargetCalculator()
      {};

      //! Just delegates to the base class function
      //! TargetCalculator::compute_default_target_matrices
    virtual void compute_target_matrices(PatchData& pd, MsqError& err)
    { compute_default_target_matrices(pd, err); MSQ_CHKERR(err); }

  protected:
    
  private:
    
  };


} //namespace


#endif // DefaultTargetCalculator_hpp
