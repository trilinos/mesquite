// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 14-Nov-02 at 17:33:19 by Thomas Leurent
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
#include "PlanarDomain.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "MesquiteError.hpp"
#include "MeshSet.hpp"

// algorythms
#include "MeanRatioQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MsqMessage.hpp"


using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "main"
int main()
{
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  MsqError err;
    //create geometry: plane z=0, normal (0,0,1)
  Vector3D pnt(0,0,0);
  Vector3D s_norm(0,0,1);
  Mesquite::PlanarDomain msq_geom(s_norm, pnt, mesh);
     
  mesh->read_vtk("../../meshFiles/2D/VTK/hybrid_3quad_1tri.vtk", err);
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.set_domain_constraint(&msq_geom);
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = MeanRatioQualityMetric::create_new();
    //   mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    //   mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
    //mean_ratio->set_averaging_method(QualityMetric::SUM, err);
    //MSQ_CHKERR(err);
  
    // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 2, err);
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  
    // creates the steepest descent, feas newt optimization procedures
    //ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  
  QualityAssessor qa=QualityAssessor(mean_ratio,QualityAssessor::ALL_MEASURES);
  
    // **************Set termination criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
    //_inner.add_criterion_type_with_double(
    // TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, OF_value, err);
    //tc_inner.add_criterion_type_with_double(
    //TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE, OF_value, err);
  TerminationCriterion tc_outer;
    //tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

    // sets a culling method on the first QualityImprover
    //This is an old command that still needs to be there.  It has
    //nothing to do with 'culling methods' described in TerminationCriterion.
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  queue1.add_quality_assessor(&qa,err); MSQ_CHKERR(err);
    // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  queue1.add_quality_assessor(&qa,err); MSQ_CHKERR(err);
  mesh->write_vtk("original_mesh",err); MSQ_CHKERR(err);
  
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  mesh->write_vtk("smoothed_mesh",err); MSQ_CHKERR(err);
  
  PRINT_TIMING_DIAGNOSTICS();
}
 
