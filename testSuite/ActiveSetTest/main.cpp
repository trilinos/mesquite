// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 17-Jul-02 at 13:53:05 by Thomas Leurent
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

#ifdef MESQUITE_USES_TSTT
#include "TSTT_Base.h"
#else
#include "TSTT_C.h"
#endif

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"

#include "CubitMesh.h" //  for writeShowMeMesh()

// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "MinTemplate.hpp"
#include "NonSmoothSteepestDescent.hpp"

using namespace Mesquite;


int main()
{     
  /* Reads a Mesh file */
  Mesh_Handle mesh;

  char file_name[128];
  strcpy(file_name, "../MeshFiles/CUBIT/2D/bad_circle");
  TSTT_Mesh_loadFile(&mesh, file_name);
  
  // Mesquite error object
  MsqError err;
  
  // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  mesh_set1.set_space_dim(3);
  mesh_set1.set_element_type(Mesquite::TRIANGLE);

  // Creates an intruction queue
  InstructionQueue queue1;

  // Creates a condition number quality metric 
  ShapeQualityMetric* cond_no = ConditionNumberQualityMetric::create_new();

  // Build an objective function with the quality metric
  MinTemplate* obj_func_min = new MinTemplate(cond_no);
  
  // Create the NonSmooth Steepest Descent procedures
  NonSmoothSteepestDescent *maxmin_method = new NonSmoothSteepestDescent( obj_func_min );
  
  // Set a culling method on the first QualityImprover
  maxmin_method->add_culling_method(QualityImprover::NO_BOUNDARY_VTX);

  // Set the max min method to be the master quality improver
  queue1.set_master_quality_improver(maxmin_method, err); MSQ_CHKERR(err);

  // write out the original mesh
  writeShowMeMesh("original_mesh", mesh);

  // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);

  // write out the smoothed mesh
  writeShowMeMesh("smoothed_mesh", mesh);

}
