// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
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
/*! \file sphere_example.cpp

This is a similar example to main.cpp for laplacian smoothing, except
the mesh is much larger and the stopping criterion is set to 500
iterations.

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
#include "MesquiteUtilities.hpp" //  for writeShowMeMesh()
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "StoppingCriterion.hpp"
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
  char file_name[128];
  /* Reads a TSTT Mesh file */
  TSTT::Mesh_Handle mesh;
  TSTT::MeshError tstt_err;
  TSTT::Mesh_Create(&mesh, &tstt_err);
  strcpy(file_name, "../../meshFiles/2D/VTK/Mesquite_geo_10242.vtk");
  TSTT::Mesh_Load(mesh, file_name, &tstt_err);
  
  // Mesquite error object
  MsqError err;
  
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
 stop_qa.add_quality_assessment(lapl_met,QualityAssessor::AVERAGE,err);
 
   //**************Set stopping criterion****************
 StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,500);
 lapl1->set_stopping_criterion(&sc2);
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
  
  writeVtkMesh("smoothed_mesh", mesh, err); MSQ_CHKERR(err);

}
