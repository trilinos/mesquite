// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD:  7-Nov-02 at 15:39:47 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifdef USE_STD_INCLUDES
#include <iostream>
#else
#include <iostream.h>
#endif

#ifdef USE_C_PREFIX_INCLUDES
#include <cstdlib>
#else
#include <stdlib.h>
#endif


#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MesquiteUtilities.hpp" //  for writeShowMeMesh()
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorithms
#include "ConditionNumberQualityMetric.hpp"
#include "MaxTemplate.hpp"
#include "NonSmoothSteepestDescent.hpp"

#include "MeshImpl.hpp"
using namespace Mesquite;


int main()
{     
    /* Reads a Mesh file */
  const char *file_name = 
//      "../../meshFiles/2D/VTK/equil_tri2.vtk";
//      "../../meshFiles/2D/VTK/bad_circle_tri_rhr.vtk";
//      "../../meshFiles/2D/VTK/tri_20258.vtk";
//      "../../meshFiles/3D/VTK/tet_1.vtk";
//      "../../meshFiles/3D/VTK/cube_tet_2.vtk";
     "../../meshFiles/3D/VTK/tire.vtk";
  printf("Loading mesh set 1\n");
  MsqError err;
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk(file_name, err);
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
    //  printf("Creating mesh set 1\n");
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
    // Creates an intruction queue
    //  printf("Creating instruction queue\n");
  InstructionQueue queue1;

  // Creates a condition number quality metric 
  //  printf("Creating quality metric\n");
  ShapeQualityMetric* cond_no = ConditionNumberQualityMetric::create_new();

  // Build an objective function with the quality metric
  //  printf("min template\n");
  MaxTemplate obj_func_min(cond_no);
  
  // Create the NonSmooth Steepest Descent procedures
  //  printf("creating optimizer\n");
  NonSmoothSteepestDescent minmax_method( &obj_func_min );
   
  minmax_method.set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err, 1);

  // Set a culling method on the first QualityImprover
  minmax_method.add_culling_method(PatchData::NO_BOUNDARY_VTX);

  // Set a termination criterion
  TerminationCriterion tc2;
  tc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  minmax_method.set_outer_termination_criterion(&tc2);
  // Set up the quality assessor
  //  printf("Setting up the quality assessor\n");
  QualityAssessor quality_assessor=QualityAssessor(cond_no,QualityAssessor::MAXIMUM);
  quality_assessor.add_quality_assessment(cond_no,QualityAssessor::MINIMUM, err); 
                   MSQ_CHKERR(err);
  quality_assessor.add_quality_assessment(cond_no,QualityAssessor::AVERAGE, err);
                   MSQ_CHKERR(err);

  // assess the quality of the initial mesh
  queue1.add_quality_assessor(&quality_assessor, err); MSQ_CHKERR(err);

  // Set the max min method to be the master quality improver
  queue1.set_master_quality_improver(&minmax_method, err); MSQ_CHKERR(err);

  // assess the quality of the final mesh
  queue1.add_quality_assessor(&quality_assessor, err); MSQ_CHKERR(err);

  // write out the original mesh
  //  printf("Writing out the original mesh\n");
  mesh->write_vtk("original_mesh", err); MSQ_CHKERR(err);

  // launches optimization on mesh_set1
  //  printf("Running the instruction queue\n");
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);

  // write out the smoothed mesh
  //  printf("Writing out the final mesh\n");
  mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
  delete cond_no;
}
