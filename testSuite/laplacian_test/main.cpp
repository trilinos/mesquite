// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD:  3-Dec-02 at 11:44:35 by Thomas Leurent
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
#include "MesquiteError.hpp"
#include "MeshImpl.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "LaplacianSmoother.hpp"
#include "EdgeLengthQualityMetric.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "main"
int main()
{     
    /* Read a VTK Mesh file */
  MsqError err;
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk("../../meshFiles/2D/VTK/square_quad_2.vtk", err);
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ShapeQualityMetric* shape_metric = ConditionNumberQualityMetric::create_new();
  SmoothnessQualityMetric* lapl_met = EdgeLengthQualityMetric::create_new();
  lapl_met->set_averaging_method(QualityMetric::RMS,err);
  
    // creates the laplacian smoother  procedures
  LaplacianSmoother* lapl1 = new LaplacianSmoother(err);
  QualityAssessor stop_qa=QualityAssessor(shape_metric,QualityAssessor::MAXIMUM);
  stop_qa.add_quality_assessment(lapl_met,QualityAssessor::ALL_MEASURES,err);
  
    //**************Set stopping criterion****************
  TerminationCriterion sc2;
  sc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,10,err);
  lapl1->set_outer_termination_criterion(&sc2);
    // sets a culling method on the first QualityImprover
  lapl1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  
    // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&stop_qa,err); MSQ_CHKERR(err);
  queue1.set_master_quality_improver(lapl1, err); MSQ_CHKERR(err);
  queue1.add_quality_assessor(&stop_qa,err); MSQ_CHKERR(err);
    // adds 1 passes of pass2 to mesh_set1
    //  mesh_set1.add_quality_pass(pass2);
  
    //writeVtkMesh("original_mesh", mesh, err); MSQ_CHKERR(err);
  
    // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
  mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
}
