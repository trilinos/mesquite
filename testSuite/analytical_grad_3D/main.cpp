// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 22-May-03 at 09:08:52 by Michael Brewer
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
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "MeanRatioQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"

using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "main"
int main()
{
  MsqError err;
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk("../../meshFiles/3D/VTK/hexes_4by2by2.vtk", err);
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
    // dbg
  std::cout << " TSTT mesh handle: " << mesh << std::endl;
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = MeanRatioQualityMetric::create_new();
  ShapeQualityMetric* cond_num = ConditionNumberQualityMetric::create_new();
  mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
  mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
//  mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  
    // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(mean_ratio);
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 2, err);
//   obj_func->set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
   // creates the steepest descent optimization procedures
  SteepestDescent* pass1 = new SteepestDescent( obj_func );
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  pass1->set_maximum_iteration(6);
  
  QualityAssessor stop_qa=QualityAssessor(mean_ratio,QualityAssessor::MAXIMUM);
  stop_qa.add_quality_assessment(cond_num, QualityAssessor::MAXIMUM,err);
  
  
   //**************Set stopping criterion****************
// StoppingCriterion sc1(&stop_qa,1.0,1.8);
    //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,1);
  TerminationCriterion tc2;
  tc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
// CompositeAndStoppingCriterion sc(&sc1,&sc2);
  pass1->set_inner_termination_criterion(&tc2);
 // sets a culling method on the first QualityImprover
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);

  // adds 1 pass of pass1 to mesh_set1
//  queue1.add_preconditioner(pass1, err); MSQ_CHKERR(err);
  queue1.add_quality_assessor(&stop_qa,err);
  queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  queue1.add_quality_assessor(&stop_qa,err);
  // adds 1 passes of pass2 to mesh_set1
//  mesh_set1.add_quality_pass(pass2);

  mesh->write_vtk("original_mesh", err); MSQ_CHKERR(err);
  
    // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
  mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
}
