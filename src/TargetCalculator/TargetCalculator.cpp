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
  pd.allocate_corner_matrices(err); MSQ_CHKERR(err);
    
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
  pd.allocate_corner_matrices(err); MSQ_CHKERR(err);
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




/*! This function computes the \f$ \lambda \f$ coefficient when it is Mesh-Based, i.e. depends on
  the whole mesh (e.g. an average of the value on the whole mesh).
  If called whereas \f$ \lambda \f$ is set to an element-based type, this function returns 1.0, so
  that is can be called safely in any context.
  
  See also the TargetCalculator::compute_Lk function for use when the \f$ \lambda \f$ coefficient
  is element-based, i.e depends only on the geometry of one element.
 */
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_L" 
double TargetCalculator::compute_L(enum lambda_type l_type, MsqError &err)
{
    // Returns 1.0 if a local (element-based) lambda was needed. 
    switch (l_type) {
    case L11:
    case L12:
    case L13:
    case L21:
    case L22:
      break;
    case L31:
    case L32:
    case L41:
      return 1.0;
    }

    // Creates a Global Patch.
    PatchDataParameters pd_params;
    PatchData ref_pd;
    pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err); MSQ_CHKERR(err);
    refMesh->get_next_patch(ref_pd, pd_params, err ); MSQ_CHKERR(err);

    // Computes lambda 
    MsqMeshEntity* elems = ref_pd.get_element_array(err); MSQ_CHKERR(err);
    double lambda=0. ;
    Matrix3D W_tri, W_quad, W_tet, W_hex;
    initialize_default_target_matrices(W_tri, W_quad, W_tet, W_hex);
//     double det_W_tri=det(W_tri);
//     double det_W_quad=det(W_quad);
    double det_W_tet=det(W_tet);
    double det_W_hex=det(W_hex);
    size_t nb_corners=0;
    Matrix3D corners[MSQ_MAX_NUM_VERT_PER_ENT];

    switch (l_type) {
    case L22:
      for (size_t i=0; i<ref_pd.num_elements(); ++i) {
        switch( elems[i].get_element_type() ) {
        case TETRAHEDRON:
          nb_corners += 4;
          elems[i].compute_corner_matrices(ref_pd, corners, 4, err); MSQ_CHKERR(err);
          for (int c=0; c<4; ++c)
            lambda += det(corners[c]) / det_W_tet;
          break;
        case HEXAHEDRON:
          nb_corners += 8;
          elems[i].compute_corner_matrices(ref_pd, corners, 8, err); MSQ_CHKERR(err);
          for (int c=0; c<8; ++c)
            lambda += det(corners[c]) / det_W_hex;
          break;
        default:
          err.set_msg("L22 not implemented for element type.");
          return 0;
        }
      }
      return pow(lambda/nb_corners, 1/3);
      
    default:
      err.set_msg("Lambda type not implemented");
      return 0;
    }

  }
    

/*! This function computes the \f$ \lambda \f$ coefficient when it is element-based, , i.e. it
  depends only on the geometry of one element.
  If called whereas \f$ \lambda \f$ is set to a mesh-based type, this function returns 1.0, so
  that is can be called safely in any context.
  
  See also the TargetCalculator::compute_L function for use when the \f$ \lambda \f$ coefficient
  is mesh-based, i.e. depends on  the whole mesh.
 */
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_Lk" 
void TargetCalculator::compute_Lk(enum lambda_type l_type, PatchData &ref_pd,
                                   size_t elem_ind, double L_k[], int num, MsqError &err)
  {
    // Returns 1.0 if a local (element-based) lambda was needed. 
    switch (l_type) {
    case L11:
    case L12:
    case L13:
    case L21:
    case L22:
      for (int i=0; i<num; ++i)
        L_k[i] = 1.0;
      return;
    case L31:
    case L32:
    case L41:
      err.set_msg("Lambda type not implemented yet");
      return;
    }

  }

  
/*! This function computes the \f$ D \f$ Matrix when it is Mesh-Based, i.e. depends on
  the whole mesh (e.g. an average of the value on the whole mesh).
  If called whereas \f$ D \f$ is set to an element-based type, this function returns
  a 3*3 identity matrix, so that is can be called safely in any context.
  
  See also the TargetCalculator::compute_Dk function for use when the \f$ D \f$ diagonal matrix
  is element-based, i.e depends only on the geometry of one element.
 */
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_D" 
 Matrix3D TargetCalculator::compute_D(enum D_type type, MsqError &err)
  {
    // Returns 1.0 if a local (element-based) lambda was needed. 
    switch (type) {
    case D11:
    case D21:
    case D31:
      break;
    case D41:
    case D42:
    case D43:
    case D51:
    case D52:
    case D53:
      Matrix3D id;
      id[0][0]=1; id[1][1]=1; id[2][2]=1;
      return id;
    }

    // Creates a Global Patch.
    PatchDataParameters pd_params;
    PatchData ref_pd;
    pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err); MSQ_CHKERR(err);
    refMesh->get_next_patch(ref_pd, pd_params, err ); MSQ_CHKERR(err);

    // Compute mesh-based D 
    switch (type) {
    default:
      err.set_msg("D type not implemented");
      Matrix3D zero;
      return zero;
    }

  }
    

/*! This function computes the \f$ D \f$ diagonal matrix when it is element-based, , i.e. it
  depends only on the geometry of one element.
  If called whereas \f$ D \f$ is set to a mesh-based type, this function returns a 3*3 identity
  matrix, so that is can be called safely in any context.
  
  See also the TargetCalculator::compute_D function for use when the \f$ D \f$ coefficient
  is mesh-based, i.e. depends on  the whole mesh.
 */
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_Dk" 
 void TargetCalculator::compute_Dk(enum D_type type, PatchData &ref_pd, size_t elem_ind,
                                           Matrix3D D_k[], int num, MsqError &err)
  {
    // Returns 1.0 if a global (mesh-based) D was needed. 
    switch (type) {
    case D11:
    case D21:
    case D31:
      for (int i=0; i<num; ++i) {
        D_k[i].zero();
        D_k[i][0][0]=1; D_k[i][1][1]=1; D_k[i][2][2]=1;
      }
      return;
    case D41:
    case D42:
    case D43:
    case D51:
    case D52:
    case D53:
      break;
    }

    MsqMeshEntity* elems = ref_pd.get_element_array(err); MSQ_CHKERR(err);
    Matrix3D corners[MSQ_MAX_NUM_VERT_PER_ENT];
    size_t nve = elems[elem_ind].vertex_count();
    
    switch(type) {
    case D31:
      elems[elem_ind].compute_corner_matrices(ref_pd, corners, nve, err); MSQ_CHKERR(err);
      for (size_t c=0; c<nve; ++c) {
        D_k[c].zero();
        D_k[c][0][0] = corners[c].column_length(0);
        D_k[c][1][1] = corners[c].column_length(1);
        D_k[c][2][2] = corners[c].column_length(2);
      }
    default:
      err.set_msg("D type not implemented.");
    }
    
  }

  
/*! This function computes the \f$ R \f$ Matrix when it is Mesh-Based, i.e. depends on
  the whole mesh (e.g. an average of the value on the whole mesh).
  If called whereas \f$ R \f$ is set to an element-based type, this function returns
  a 3*3 identity matrix, so that is can be called safely in any context.
  
  See also the TargetCalculator::compute_Rk function for use when the \f$ R \f$ diagonal matrix
  is element-based, i.e depends only on the geometry of one element.
 */
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_R" 
 Matrix3D TargetCalculator::compute_R(enum R_type type, MsqError &err)
  {
    // Returns 1.0 if a local (element-based) lambda was needed. 
    switch (type) {
    case R11:
    case R21:
    case R31:
    case R32:
    case R33:
      break;
    case R00:
    case R41:
    case R42:
    case R43:
    case R44:
      Matrix3D id;
      id[0][0]=1; id[1][1]=1; id[2][2]=1;
      return id;
    }

    // Creates a Global Patch.
    PatchDataParameters pd_params;
    PatchData ref_pd;
    pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err); MSQ_CHKERR(err);
    refMesh->get_next_patch(ref_pd, pd_params, err ); MSQ_CHKERR(err);

    // Compute mesh-based R 
    switch (type) {
    default:
      err.set_msg("R type not implemented");
      Matrix3D zero;
      return zero;
    }

  }
    

/*! This function computes the \f$ R \f$ rotation matrix when it is element-based, , i.e. it
  depends only on the geometry of one element.
  If called whereas \f$ R \f$ is set to a mesh-based type, this function returns a 3*3 identity
  matrix, so that is can be called safely in any context.
  
  See also the TargetCalculator::compute_R function for use when the \f$ R \f$ coefficient
  is mesh-based, i.e. depends on  the whole mesh.
 */
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_Rk" 
 void TargetCalculator::compute_Rk(enum R_type type, PatchData &ref_pd, size_t elem_ind,
                                           Matrix3D R_k[], int num, MsqError &err)
  {
    // Returns 1.0 if a global (mesh-based) R was needed. 
    switch (type) {
    case R00:
    case R11:
    case R21:
    case R31:
    case R32:
    case R33:
      for (int i=0; i<num; ++i) {
        R_k[i].zero();
        R_k[i][0][0]=1; R_k[i][1][1]=1; R_k[i][2][2]=1;
      }
      return;
    case R41:
    case R42:
    case R43:
    case R44:
      break;
    }

    switch(type) {
    case W00:
      
    default:
      err.set_msg("R type not implemented.");
    }
    
  }


/*! This function computes the \f$ W \f$  matrix when it is element-based, , i.e. it
  depends only on the geometry of one element.
  If called whereas \f$ R \f$ is set to a mesh-based type, this function returns a 3*3 identity
  matrix, so that is can be called safely in any context.
*/
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_Wk" 
 void TargetCalculator::compute_Wk(enum W_type type, PatchData &ref_pd, size_t elem_ind,
                                           Matrix3D W_k[], int num, MsqError &err)
  {
    // Set matrices to identity if a global (mesh-based) W was needed. 
    switch (type) {
    case W00:
    case W11:
    case W21:
    case W31:
    case W41:
    case W42:
      break;
    }
    
    MsqMeshEntity* elems = ref_pd.get_element_array(err); MSQ_CHKERR(err);
    size_t nve = elems[elem_ind].vertex_count();

    switch(type) {
    case W00:
      elems[elem_ind].compute_corner_matrices(ref_pd, W_k, nve, err); MSQ_CHKERR(err);
    default:
      err.set_msg("W type not implemented.");
    }
    
  }



/*! The type of targets computed by this function is selected by the constructor of
    the base classes. */
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
    compute_Lk(mLambda, ref_pd, i, Lk, nve, err); MSQ_CHKERR(err);
    compute_Dk(mD, ref_pd, i, Dk, nve, err); MSQ_CHKERR(err);
    compute_Rk(mR, ref_pd, i, Rk, nve, err); MSQ_CHKERR(err);
    for (int i=0; i<nve; ++i) {
      tag->target_matrix(i) = L * D * R * Lk[i] * Rk[i] * Dk[i] * W[i];
    }
  }
    
  FUNCTION_TIMER_END();
}

