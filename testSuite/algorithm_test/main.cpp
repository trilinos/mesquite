// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD:  2-Apr-03 at 17:23:53 by Thomas Leurent
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

// algorythms
#include "MeanRatioQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "main"
int main()
{     
  MsqError err;
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk("../../meshFiles/3D/VTK/tire.vtk", err);
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
    // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean = MeanRatioQualityMetric::create_new();
//   mean->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
//   mean->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
  mean->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
  mean->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
  
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean, 1, err);
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
    //obj_func->set_hessian_type(ObjectiveFunction::ANALYTICAL_HESSIAN);
  
  // creates the optimization procedures
//   ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func );

  //perform optimization globally
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
  
  QualityAssessor mean_qa=QualityAssessor(mean,QualityAssessor::AVERAGE);

    //**************Set termination criterion****************

  //perform 1 pass of the outer loop (this line isn't essential as it is
  //the default behavior).
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES, 1, err);
  pass1->set_outer_termination_criterion(&tc_outer);
  
  //perform the inner loop until a certain objective function value is
  //reached.  The exact value needs to be determined (about 18095).
  //As a safety, also stop if the time exceeds 10 minutes (600 seconds).
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, 13975, err);
//  tc_inner.add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, 13964.93818, err);
  tc_inner.add_criterion_type_with_double(TerminationCriterion::CPU_TIME, 1800, err);
  
  pass1->set_inner_termination_criterion(&tc_inner);
  
    // sets a culling method on the first QualityImprover
  //This should probably be removed
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  //used for cg to get some info
//  pass1->set_debugging_level(2);
  
  // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&mean_qa,err);
  queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  queue1.add_quality_assessor(&mean_qa,err);
  mesh->write_vtk("original_mesh", err); MSQ_CHKERR(err);
  
    // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
  mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
  PRINT_TIMING_DIAGNOSTICS();
}
