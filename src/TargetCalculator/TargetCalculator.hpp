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
    
    //! virtual destructor ensures use of polymorphism during destruction
    virtual ~TargetCalculator()
      {};

      //! Compute the default "isotropic" target matrices that are often used in the computation
      //! of reference-based target matrices.
      //! The resulting corner matrices are stored in tags on the elements of the PatchData.
    void compute_default_target_matrices(PatchData &pd, MsqError &err);
    
    /*! \brief This function must provide the corner matrices for all elements on the Patch.

         Useful functionality includes: MsqMeshEntity::set_tag, MsqTag::target_matrix,
         MsqTag::scalar . 
    */
    virtual void compute_target_matrices(PatchData& pd, MsqError& err)=0;

  protected:
    
  private:
    
  };


#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_default_target_matrices" 
  inline void TargetCalculator::compute_default_target_matrices(PatchData &pd,
                                                      MsqError &err)
  {
    FUNCTION_TIMER_START(__FUNC__);
    
    // set on each element in the patch a tag containing an array of corner matrices
    // (the size of the array is adequate for each element, e.g. 4 for a quad).
    pd.allocate_corner_matrices(err); MSQ_CHKERR(err);
    
    MsqMeshEntity* elems=pd.get_element_array(err);
    size_t num_elements=pd.num_elements();

    const double v_tri[] = {1, 0.5, 0, 0, MSQ_SQRT_THREE/2, 0, 0, 0, 0};
    Matrix3D tmp_tri(v_tri);

    const double v_quad[] = {1, 0, 0, 0, 1, 0, 0, 0, 0};
    Matrix3D tmp_quad(v_quad);
    
    const double v_tet[] = {1, 0.5, 0.5, 0, MSQ_SQRT_THREE/2, MSQ_SQRT_THREE/6, 0, 0, MSQ_SQRT_TWO/MSQ_SQRT_THREE};
    Matrix3D tmp_tet(v_tet);

    const double v_hex[] = {1, 0, 0,  0, 1, 0,  0, 0, 1};
    Matrix3D tmp_hex(v_hex);

    // set the corner matrices to the correct value for each tag.
    for (size_t i=0; i<num_elements; ++i) {

      MsqTag* tag = elems[i].get_tag(); 
      
      EntityTopology type = elems[i].get_element_type();
      switch (type)
        {
        case TRIANGLE:
          tag->target_matrix(0) = tmp_tri; 
          tag->target_matrix(1) = tmp_tri; 
          tag->target_matrix(2) = tmp_tri; 
          break;
        case QUADRILATERAL:
          tag->target_matrix(0) = tmp_quad; 
          tag->target_matrix(1) = tmp_quad; 
          tag->target_matrix(2) = tmp_quad; 
          tag->target_matrix(3) = tmp_quad; 
          break;
        case TETRAHEDRON:
          tag->target_matrix(0) = tmp_tet; 
          tag->target_matrix(1) = tmp_tet; 
          tag->target_matrix(2) = tmp_tet; 
          tag->target_matrix(3) = tmp_tet; 
          break;
        case HEXAHEDRON:
          tag->target_matrix(0) = tmp_hex; 
          tag->target_matrix(1) = tmp_hex; 
          tag->target_matrix(2) = tmp_hex; 
          tag->target_matrix(3) = tmp_hex; 
          tag->target_matrix(4) = tmp_hex; 
          tag->target_matrix(5) = tmp_hex; 
          tag->target_matrix(6) = tmp_hex; 
          tag->target_matrix(7) = tmp_hex; 
          break;
        default:
          err.set_msg("Type not implemented.");
          return;
        } //end switch
    } // end loop
    FUNCTION_TIMER_END();   
  }

  
  
} //namespace


#endif // TargetCalculator_hpp
