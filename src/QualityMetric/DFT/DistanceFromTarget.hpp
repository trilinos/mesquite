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

/*! \file DistanceFromTarget.hpp

Header file for the Mesquite::DistanceFromTarget class

  \author Thomas Leurent
  \date   2004-09-29
 */


#ifndef DistanceFromTarget_hpp
#define DistanceFromTarget_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MsqMessage.hpp"
#include "QualityMetric.hpp"
#include "MsqMeshEntity.hpp"

namespace Mesquite
{
  
  /*! \class DistanceFromTarget
    \brief Base class for the computation of the distance from target between
           the target matrices W and the actual corner matrices A. 
  */
  class DistanceFromTarget : public QualityMetric
  {
  public:
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~DistanceFromTarget()
       {};

  protected:
      //! For a given element, compute each corner matrix A, and given a target
      //! corner matrix W, returns \f$ T=AW^{-1} \f$ for each corner.    
    void compute_T_matrices(MsqMeshEntity &elem, PatchData& pd,
                            Matrix3D T[], size_t num_T, MsqError &err);
 
  private:
    
  };

  inline void DistanceFromTarget::compute_T_matrices(MsqMeshEntity &elem, PatchData& pd,
                                 Matrix3D T[], size_t num_T, MsqError &err)
  {    
      // Gets the element corner matrices.
    elem.compute_corner_matrices(pd, T, num_T, err);

//     for (size_t i=0; i<num_T; ++i)
//       std::cout << "A["<<i<<"]:\n" << T[i] << std::endl;

    size_t num_targets=0;
    TargetMatrix* W = elem.get_target_matrices(num_targets, err); MSQ_CHKERR(err);
    assert( num_T == num_targets );
    
//     for (size_t i=0; i<num_T; ++i)
//       std::cout << "W["<<i<<"]:\n" << W[i] << std::endl;

    for (size_t i=0; i<num_T; ++i)
      timesInvA(T[i], W[i]);
  }

  
  
} //namespace


#endif // DistanceFromTarget_hpp
