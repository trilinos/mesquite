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

/*! \file TargetCalculator.cpp

\brief The Mesquite::TargetCalculator class is the base class. Concrete classes are 
instantiated by the user, and often implemented by the user to give 
mesquite a measure of the perfect mesh. 

\author Thomas Leurent
\date   2004-04-30
*/


#include "TargetCalculator.hpp"
#include "PatchDataUser.hpp"
#include "MeshSet.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_target_matrices_and_check_det" 
void TargetCalculator::compute_target_matrices_and_check_det(PatchData &pd, MsqError &err)
{
  FUNCTION_TIMER_START(__FUNC__);

  // Compute the target matrices
  compute_target_matrices(pd, err); MSQ_CHKERR(err);

  //checks that the determinant of each target matrix is positive.
  MsqMeshEntity* elems=pd.get_element_array(err);
  size_t num_elements=pd.num_elements();
  for (size_t i=0; i<num_elements; ++i) {
    MsqTag* tag = elems[i].get_tag();
    size_t num_corners = elems[i].vertex_count();
    for (size_t j=0; j<num_corners; ++j) {    
      if ( det(tag->target_matrix(j)) <= 0 ) {
        err.set_msg("A Target matrix has a non-positive determinant. Please review your target calculator.");
        FUNCTION_TIMER_END();
        return;
      }
    }
  }
    
  FUNCTION_TIMER_END();
}

  
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_default_target_matrices" 
void TargetCalculator::compute_default_target_matrices(PatchData &pd,
                                                       MsqError &err)
{
  FUNCTION_TIMER_START(__FUNC__);
    
  // set on each element in the patch a tag containing an array of corner matrices
  // (the size of the array is adequate for each element, e.g. 4 for a quad).
  pd.allocate_target_matrices(err); MSQ_CHKERR(err);
    
  MsqMeshEntity* elems=pd.get_element_array(err);
  size_t num_elements=pd.num_elements();

  Matrix3D tmp_tri, tmp_quad, tmp_tet, tmp_hex;
  initialize_default_target_matrices(tmp_tri, tmp_quad, tmp_tet, tmp_hex);
  
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

  
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_reference_corner_matrices" 
void TargetCalculator::compute_reference_corner_matrices(PatchData &pd,
                                                         MsqError &err)
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
    
  // Compute corner matrices for ref_pd and store in pd tags.
  MsqMeshEntity* elems = pd.get_element_array(err);
  MsqMeshEntity* elems_ref = ref_pd.get_element_array(err);
  TargetMatrix A[MSQ_MAX_NUM_VERT_PER_ENT];
  pd.allocate_target_matrices(err); MSQ_CHKERR(err);
  for (size_t i=0; i<num_elements; ++i) {
    MsqTag* tag = elems[i].get_tag();
    int nve = elems[i].vertex_count();
    assert( nve = elems_ref[i].vertex_count() );
    elems_ref[i].compute_corner_matrices(ref_pd, A, nve, err);
    for (int i=0; i<nve; ++i) {
      tag->target_matrix(i) = A[i];
    }
  }
    
  FUNCTION_TIMER_END();
}


#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_guide_matrix" 
 void TargetCalculator::compute_guide_matrices(enum guide_type type, PatchData &ref_pd, size_t elem_ind,
                                           Matrix3D W_k[], int num, MsqError &err)
  {

    
    MsqMeshEntity* elems = ref_pd.get_element_array(err); MSQ_CHKERR(err);
    size_t nve = elems[elem_ind].vertex_count();

    switch(type) {
    case Ad:
      {
        Matrix3D tmp_tri, tmp_quad, tmp_tet, tmp_hex;
        initialize_default_target_matrices(tmp_tri, tmp_quad, tmp_tet, tmp_hex);
        EntityTopology elem_type = elems[elem_ind].get_element_type();
        switch (elem_type) {
        case TRIANGLE:
          for (int i=0; i<3; ++i) W_k[i] = tmp_tri; 
          break;
        case QUADRILATERAL:
          for (int i=0; i<4; ++i) W_k[i] = tmp_quad; 
          break;
        case TETRAHEDRON:
          for (int i=0; i<4; ++i) W_k[i] = tmp_tet; 
          break;
        case HEXAHEDRON:
          for (int i=0; i<8; ++i) W_k[i] = tmp_hex; 
          break;
        default:
          err.set_msg("Element type not implemented.");
          return;
        }
      }
      return;
    case A0:
      elems[elem_ind].compute_corner_matrices(ref_pd, W_k, nve, err); MSQ_CHKERR(err);
      return;
    case Ar:
      elems[elem_ind].compute_corner_matrices(ref_pd, W_k, nve, err); MSQ_CHKERR(err);
      return;
    default:
      err.set_msg("'A' guide type not implemented.");
    }
    
  }

