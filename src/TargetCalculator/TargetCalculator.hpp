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

/*! \file TargetCalculator.hpp

\brief The Mesquite::TargetCalculator class is the base class. Concrete classes are 
 instantiated by the user, and often implemented by the user to give 
 mesquite a measure of the perfect mesh. 

  \author Thomas Leurent
  \date   2004-09-31
 */


#ifndef TargetCalculator_hpp
#define TargetCalculator_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MsqTimer.hpp"
#include "MsqMessage.hpp"
#include "PatchData.hpp"

namespace Mesquite
{
  class PatchDataParameters;

  
  /*! \class TargetCalculator
    \brief Base class that provides the interface for computing the target corner matrices
    used in the context based smoothing.

    To implement a concrete TargetCalculator, one inherits from this class and then 
    overrides the compute_target_matrices function itself to provide corner matrices for all
    elements in the patch.

    Note that an implementation is provided in TargetCalculator::compute_default_target_matrices
    for default target matrices often used in the computation of more complex,
    reference-based target matrices.

    The target calculator is set on the QualityImprover. At runtime, it associates with
    each mesh corner (i.e. a corner of an element. Most of the time, one vertex
    corresponds to many corners.) with a TargetMatrix through an MsqTag mechanism 
    available in the MsqMeshEntity class.
   */
  class TargetCalculator 
  {
  public:

    TargetCalculator() : refMesh(0), originator(0) { }
    
    //! virtual destructor ensures use of polymorphism during destruction
      virtual ~TargetCalculator()
    {};

      //! Compute the default "isotropic" target matrices that are often used in the computation
      //! of reference-based target matrices.
      //! The resulting corner matrices are stored in tags on the elements of the PatchData.
    void compute_default_target_matrices(PatchData &pd, MsqError &err);


      //! Compute the corner matrices for the reference mesh refMesh.
      //! The refMesh data member is set by the constructors of a concrete TargetCalculator
      //! that requires a reference mesh.
    void compute_reference_corner_matrices(PatchData &pd, MsqError &err);

    //! This function wraps compute_target_matrices and checks that the determinant of each target
    //! is positive.
    void compute_target_matrices_and_check_det(PatchData& pd, MsqError& err);
    
    /*! \brief This function must provide the corner matrices for all elements on the Patch.

         Useful functionality includes: MsqMeshEntity::set_tag, MsqTag::target_matrix,
         MsqTag::scalar .
    */
    virtual void compute_target_matrices(PatchData& pd, MsqError& err)=0;

    void set_originator(PatchDataParameters* pdm, MsqError &err)
    { if (originator != 0)
        err.set_msg("Each TargetCalculator can be set on one object only.");
      else originator = pdm; }

  protected:
    
  private:
    MeshSet* refMesh;
    PatchDataParameters* originator; //! This is the object the TargetCalculator is attached to.
  };

    
} //namespace


#endif // TargetCalculator_hpp
