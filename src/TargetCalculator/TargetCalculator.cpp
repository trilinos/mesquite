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
#include "MsqTimer.hpp"
#include "TargetMatrix.hpp"

using namespace Mesquite;


Matrix3D TargetCalculator::get_default_target_matrix( EntityTopology type,
                                                      MsqError& err )
{
  static const double SIXTH_ROOT_OF_TWO = msq_stdc::pow(2., 1./6.);
  Matrix3D result;
  switch( type )
  {
    case TRIANGLE:
    {
#ifndef MSQ_JANUS
      const double v_tri[] = {1., 0.5, 0., 0., MSQ_SQRT_THREE/2., 0., 0., 0., 1.};
#else    
      const double v_tri[] = {1., 0.5, 0., 0., 0.8660254037844385965883021 , 0., 0., 0., 1.};
#endif    
      Matrix3D m1(v_tri);
      m1 *= MSQ_3RT_2_OVER_6RT_3;
      return m1;
    }

    case QUADRILATERAL:
    case HEXAHEDRON:
    {
      const double ident[] = { 1.0, 0.0, 0.0, 
                               0.0, 1.0, 0.0, 
                               0.0, 0.0, 1.0 };
      return Matrix3D(ident);
    }

    case TETRAHEDRON:
    {
#ifndef MSQ_JANUS
      const double v_tet[] = {1., 0.5, 0.5, 0., MSQ_SQRT_THREE/2., MSQ_SQRT_THREE/6., 0., 0., MSQ_SQRT_TWO/MSQ_SQRT_THREE};
#else 
      const double v_tet[] = {1., 0.5, 0.5, 0.,  0.8660254037844385965883021, 0.2886751345948128655294340, 0., 0., 0.8164965809277261454823815};
#endif    
      Matrix3D m3(v_tet);
      m3 *= SIXTH_ROOT_OF_TWO;
      return m3;
    }

    case PYRAMID:
    {
      const double v_pyr[] = { 1.0, 0.0, 0.5, 
                               0.0, 1.0, 0.5, 
                               0.0, 0.0, 1.0 };
      return Matrix3D(v_pyr);
    }

    case PRISM:
    {
      const double s = MSQ_3RT_2_OVER_6RT_3;
      const double h = MSQ_3RT_2_OVER_6RT_3 * MSQ_SQRT_THREE_DIV_TWO;
      const double v_wdg[] = {  s, 0.5*s, 0,  
                                0,  h,    0,
                                0,  0,    s };
      return Matrix3D(v_wdg);
    } 

    default:
      MSQ_SETERR(err)( MsqError::UNSUPPORTED_ELEMENT,
                "Unsupported element type (%d)\n",
                (int)type);
      return Matrix3D();
  }
}

void TargetCalculator::compute_target_matrices_and_check_det( PatchData &pd, 
                                                              PatchData& ref_pd,
                                                              MsqError &err)
{
  MSQ_FUNCTION_TIMER( "TargetCalculator::compute_target_matrices_and_check_det" );

  // Compute the target matrices
  compute_target_matrices(pd, ref_pd, err); MSQ_ERRRTN(err);

  //checks that the determinant of each target matrix is positive.
  MsqMeshEntity* elems=pd.get_element_array(err); MSQ_ERRRTN(err);
  size_t num_elements=pd.num_elements();
  for (size_t i=0; i<num_elements; ++i) {
    const TargetMatrix* matrices = pd.targetMatrices.get_element_corner_tags(&pd, i, err);
    MSQ_ERRRTN(err);
    size_t num_corners = elems[i].corner_count();
    
    for (size_t j=0; j<num_corners; ++j) {    
      if ( det(matrices[j]) <= 0 ) {
        MSQ_SETERR(err)("A Target matrix has a non-positive determinant. "
                        "Please review your target calculator.",
                        MsqError::INVALID_STATE);
        return;
      }
    }
  }
}

  
void TargetCalculator::compute_default_target_matrices(PatchData &pd,
                                                       MsqError &err)
{
  MSQ_FUNCTION_TIMER( "TargetCalculator::compute_default_target_matrices" );
    
  MsqMeshEntity* elems=pd.get_element_array(err); MSQ_ERRRTN(err);
  size_t num_elements=pd.num_elements();
  TargetMatrix matrices[8];
  
  // set the corner matrices to the correct value for each tag.
  for (size_t i=0; i<num_elements; ++i) {
    EntityTopology type = elems[i].get_element_type();
    Matrix3D tmp = get_default_target_matrix( type, err ); MSQ_ERRRTN(err);
    for (size_t j = 0; j < elems[i].corner_count(); ++j)
      matrices[j]= tmp;
      
    pd.targetMatrices.set_element_corner_tags( &pd, i, matrices, err ); MSQ_ERRRTN(err);
  } // end loop
}

/*  
void TargetCalculator::compute_reference_corner_matrices(PatchData &pd,
                                                         MsqError &err)
{
  MSQ_FUNCTION_TIMER( "TargetCalculator::compute_reference_corner_matrices" );

  PatchData ref_pd;
  if (refMesh)
    refMesh->get_next_patch(ref_pd, this->get_all_parameters(), err); 
  else 
    MSQ_SETERR(err)("No reference mesh", MsqError::INVALID_STATE);
  MSQ_ERRRTN(err);
  
  // Make sure topology of ref_pd and pd are equal
  size_t num_elements=pd.num_elements();
  assert( num_elements == ref_pd.num_elements() );
  size_t num_vertices=pd.num_vertices();
  assert( num_vertices == ref_pd.num_vertices() );
    
  // Compute corner matrices for ref_pd and store in pd tags.
  MsqMeshEntity* elems = pd.get_element_array(err);
  MsqMeshEntity* elems_ref = ref_pd.get_element_array(err); MSQ_ERRRTN(err);
  TargetMatrix A[MSQ_MAX_NUM_VERT_PER_ENT];
  for (size_t i=0; i<num_elements; ++i) {
    int nve = elems[i].corner_count();
    assert( nve = elems_ref[i].corner_count() );
    elems_ref[i].compute_corner_matrices(ref_pd, A, nve, err); MSQ_ERRRTN(err);
    pd.targetMatrices.set_element_corner_tags( &pd, i, A, err ); MSQ_ERRRTN(err);
  }
}
*/

 void TargetCalculator::compute_guide_matrices(enum guide_type type, PatchData &ref_pd, size_t elem_ind,
                                           Matrix3D W_k[], int num, MsqError &err)
  {

    
    MsqMeshEntity* elems = ref_pd.get_element_array(err); MSQ_ERRRTN(err);
    size_t nve = elems[elem_ind].corner_count();
 
    switch(type) {
    case Ad:
      {
        EntityTopology elem_type = elems[elem_ind].get_element_type();
        Matrix3D tmp = get_default_target_matrix( elem_type, err ); MSQ_ERRRTN(err);
        for (unsigned i = 0; i < elems[elem_ind].corner_count(); ++i)
          W_k[i] = tmp;
      }
      return;
    case A0:
      elems[elem_ind].compute_corner_matrices(ref_pd, W_k, nve, err); MSQ_ERRRTN(err);
      return;
    case Ar:
      elems[elem_ind].compute_corner_matrices(ref_pd, W_k, nve, err); MSQ_ERRRTN(err);
      return;
    default:
      MSQ_SETERR(err)("Element type not implemented.",MsqError::NOT_IMPLEMENTED);
      return;
    }

  }

  //!
  Matrix3D TargetCalculator::compute_V_3D(const Matrix3D &A, MsqError &err)
  {
    Vector3D a1(A[0][0], A[1][0], A[2][0]); 
    Vector3D a2(A[0][1], A[1][1], A[2][1]); 
    Vector3D a3(A[0][2], A[1][2], A[2][2]); 

    double a1_norm = A.column_length(0); 
    Vector3D a1_x_a2 = a1 * a2;
    double a1_x_a2_norm = a1_x_a2.length(); 
    
    Matrix3D V;
    Vector3D v1, v2, v3;
    
    double a1_norm_inv = 1./a1_norm;
    double a1_x_a2_norm_inv = 1./a1_x_a2_norm;
    if (!finite(a1_norm_inv) || !finite(a1_x_a2_norm_inv)) {
      MSQ_SETERR(err)("Numerical error", MsqError::INVALID_STATE);
    }
    
    // note : % is the dot product
    v1 = a1_norm_inv * a1;
    v2 = ((-(a1%a2) * a1) + (a1_norm*a1_norm)*a2) * a1_norm_inv * a1_x_a2_norm_inv;
    v3 = a1_x_a2_norm_inv * a1_x_a2;

    V.set_column(0, v1);
    V.set_column(1, v2);
    V.set_column(2, v3);

    return V;
  }
  
  msq_std::string TargetCalculator::get_name()
    { return "Target Calculator"; }
  
  PatchDataUser::AlgorithmType TargetCalculator::get_algorithm_type()
    { return PatchDataUser::TARGET_CALCULATOR; }
  
  double TargetCalculator::loop_over_mesh( Mesh* mesh,
                                           MeshDomain* domain,
                                           PatchData* global_patch,
                                           MsqError& err )
  {
    PatchData patch, ref_patch;
    patch.set_mesh( mesh );
    patch.set_domain( domain );
    ref_patch.set_mesh( refMesh );
    ref_patch.set_domain( refDomain );
    const bool have_ref_mesh = refMesh && (refMesh == mesh || refDomain == domain);
    
    
    if (get_patch_type() == PatchData::GLOBAL_PATCH || global_patch)
    {
      if (!global_patch)
      {
        patch.fill_global_patch( err ); MSQ_ERRZERO(err);
        global_patch = &patch;
      }
      
      PatchData* ref_patch_ptr;
      if (!have_ref_mesh)
      {
        ref_patch_ptr = global_patch;
      }
      else
      {
        ref_patch.fill_global_patch( err ); MSQ_ERRZERO(err);
        ref_patch_ptr = &ref_patch;
      }
      
      compute_target_matrices_and_check_det( *global_patch,
                                             *ref_patch_ptr,
                                             err ); MSQ_CHKERR(err);
    }
        
    else if (!have_ref_mesh)
    {
      while (patch.get_next_element_patch( err ))
      {
        MSQ_ERRZERO(err);
        compute_target_matrices_and_check_det( patch, patch, err ); MSQ_ERRZERO(err);
      }
      MSQ_CHKERR(err);
    }
    
    else
    {
      while (patch.get_next_element_patch( err ))
      {
        MSQ_ERRZERO(err);
        if (!ref_patch.get_next_element_patch(err))
        {
          MSQ_SETERR(err)("Target mesh and reference mesh are not consistant.\n",
                          MsqError::INVALID_STATE);
          return 0;
        }
        MSQ_ERRZERO(err);
        
        compute_target_matrices_and_check_det( patch, ref_patch, err ); MSQ_ERRZERO(err);
      }
      MSQ_CHKERR(err);
    }
    
    return 1.0;
  }  
  
  
