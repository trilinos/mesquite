// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 25-Apr-03 at 13:38:11 by Thomas Leurent
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
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "MeanRatioQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"

#include "MsqMessage.hpp"
using namespace Mesquite;

using std::cout;
using std::endl;

#undef __FUNC__
#define __FUNC__ "main"
int main(int argc, char* argv[])
{
  Mesquite::MsqError err;
  char file_name[128];
  double OF_value = 100;
  
  // command line arguments
  if (argc==1 || argc>3)
    cout << "meshfile name needed as argument.\n"
      "objective function value optional as 2nd argument.\n" << endl;
  else if (argc==2) {
    cout << " given 1 command line argument.\n";
    strcpy(file_name, argv[1]);
  } else if (argc==3) {
    cout << " given 2 command line arguments.\n";
    strcpy(file_name, argv[1]);
    OF_value = atof(argv[2]);
  }
  
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk(file_name, err);
  
  // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = MeanRatioQualityMetric::create_new();
  mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
  
  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 2, err);
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  
  // creates the steepest descentfeas newt optimization procedures
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  
  QualityAssessor stop_qa=QualityAssessor(mean_ratio,QualityAssessor::MAXIMUM);
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_double(
           TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, OF_value, err);
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  // sets a culling method on the first QualityImprover
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  
  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  
  mesh->write_vtk("original_mesh",err); MSQ_CHKERR(err);
  
  // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
  mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
  PRINT_TIMING_DIAGNOSTICS();
}
