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

/*! \file LRWDTargetCalculator.hpp

Header file for the Mesquite::LRWDTargetCalculator class

  \author Thomas Leurent
  \date   2004-09-31
 */


#ifndef LRWDTargetCalculator_hpp
#define LRWDTargetCalculator_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MsqMessage.hpp"
#include "TargetCalculator.hpp"

namespace Mesquite
{
  
  /*! \class LRWDTargetCalculator
    \brief This class computes uses the corner matrices of a reference mesh
    as the target matrices.
  */
  class LRWDTargetCalculator : public TargetCalculator
  {
  public:
    
    LRWDTargetCalculator(MeshSet* ref_mesh, enum lambda_type lambda, enum D_type D, enum R_type R, err)
    { refMesh = ref_mesh;
      mLambda = lambda;
      mD = D;
      mR = R;
    }

      
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~LRWDTargetCalculator()
      {};

      //! Just delegates to the base class function
      //! TargetCalculator::compute_default_target_matrices
    virtual void compute_target_matrices(PatchData& pd, MsqError& err)
    { compute_reference_corner_matrices(pd, err); MSQ_CHKERR(err); }

  protected:
    
  private:
    enum lambda_type mLambda; 
    enum D_type mD;
    enum R_type mR;
  };


#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_target_matrices" 
void TargetCalculator::compute_target_matrices(PatchData &pd, MsqError &err)
{
  FUNCTION_TIMER_START(__FUNC__);

  if (refMesh == 0) {
    err.set_msg("Reference mesh has not been set. If the target calculator uses"
                "a reference mesh, it should set it to a constructor argument.");
    return;
  }

  PatchData ref_pd;
  refMesh->get_next_patch(ref_pd, *originator, err); MSQ_CHKERR(err);
  
  // Make sure topology of ref_pd and pd are equal
  size_t num_elements=pd.num_elements();
  assert( num_elements == ref_pd.num_elements() );
  size_t num_vertices=pd.num_vertices();
  assert( num_vertices == ref_pd.num_vertices() );
    
  // .
  Matrix3D tmp_tri, tmp_quad, tmp_tet, tmp_hex;
  initialize_default_target_matrices(tmp_tri, tmp_quad, tmp_tet, tmp_hex);
  
  MsqMeshEntity* elems = pd.get_element_array(err);
  MsqMeshEntity* elems_ref = ref_pd.get_element_array(err);
  TargetMatrix W[MSQ_MAX_NUM_VERT_PER_ENT];
  pd.allocate_corner_matrices(err); MSQ_CHKERR(err);
  double L = compute_L(mLambda, err); MSQ_CHKERR(err);
  Matrix3D D = compute_D(mD, err); MSQ_CHKERR(err);
  Matrix3D R = compute_R(mR, err); MSQ_CHKERR(err);
  for (size_t i=0; i<num_elements; ++i) {
    MsqTag* tag = elems[i].get_tag();
    int nve = elems[i].vertex_count();
    assert( nve = elems_ref[i].vertex_count() );
    elems_ref[i].compute_corner_matrices(ref_pd, W, nve, err);
    double Lk[MSQ_MAX_NUM_VERT_PER_ENT];
    Matrix3D Dk[MSQ_MAX_NUM_VERT_PER_ENT];
    Matrix3D Rk[MSQ_MAX_NUM_VERT_PER_ENT];
    compute_Lk(mLambda, ref_pd, i, Lk, err); MSQ_CHKERR(err);
    compute_Dk(mD, ref_pd, i, Dk, err); MSQ_CHKERR(err);
    compute_Rk(mR, ref_pd, i, Rk, err); MSQ_CHKERR(err);
    for (int i=0; i<nve; ++i) {
      tag->target_matrix(i) = L * D * R * Lk[i] * Rk[i] * Dk[i] * W[i];
    }
  }
    
  FUNCTION_TIMER_END();
}

  

} //namespace


#endif // LRWDTargetCalculator_hpp
