// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD:  3-Dec-02 at 11:29:33 by Thomas Leurent
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
#include "TSTT_Base.h"
#include "MeshImpl.hpp"
#include "MesquiteUtilities.hpp" //  for writeShowMeMesh()
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "MeanRatioQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"

#include "MsqMessage.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "main"
int main()
{
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  Mesquite::MsqError err;
  char file_name[128];
  strcpy(file_name, "../../meshFiles/2D/VTK/square_tri_2.vtk");
  mesh->read_vtk(file_name, err);
  
  // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);

  // dbg
   std::cout << " TSTT mesh handle: " << mesh << std::endl;
  
  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = MeanRatioQualityMetric::create_new();

  // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(mean_ratio);
  LPTemplate* obj_func = new LPTemplate(mean_ratio, 2, err);
    // creates the steepest descent optimization procedures
  SteepestDescent* pass1 = new SteepestDescent( obj_func );
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  
 QualityAssessor stop_qa=QualityAssessor(mean_ratio,QualityAssessor::MAXIMUM);

   //**************Set stopping criterion****************
// StoppingCriterion sc1(&stop_qa,1.0,1.8);
   //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,1);
 TerminationCriterion sc2;
 sc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
// CompositeAndStoppingCriterion sc(&sc1,&sc2);
 pass1->set_inner_termination_criterion(&sc2);
 // sets a culling method on the first QualityImprover
 pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);

  // adds 1 pass of pass1 to mesh_set1
//  queue1.add_preconditioner(pass1, err); MSQ_CHKERR(err);
  queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  // adds 1 passes of pass2 to mesh_set1
//  mesh_set1.add_quality_pass(pass2);

    //writeVtkMesh("original_mesh", mesh, err); MSQ_CHKERR(err);
  
  // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
    //writeVtkMesh("smoothed_mesh", mesh, err); MSQ_CHKERR(err);
  PRINT_TIMING_DIAGNOSTICS();
}
