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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 23-Jul-03 at 18:09:13 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifndef MSQ_USE_OLD_IO_HEADERS
#include <iostream>
#else
#include <iostream.h>
#endif

#ifndef MSQ_USE_OLD_C_HEADERS
#include <cstdlib>
#else
#include <stdlib.h>
#endif


#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
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


int main()
{
  MsqPrintError err(msq_stdio::cout);
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk("../../meshFiles/3D/VTK/hexes_4by2by2.vtk", err);
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); 
  if (err) return 1;
  
    // dbg
  std::cout << " TSTT mesh handle: " << mesh << std::endl;
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric(err);
  if (err) return 1;
  ShapeQualityMetric* cond_num = new ConditionNumberQualityMetric;
  mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
  if (err) return 1;
  mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
//  mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  
    // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(mean_ratio);
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 2, err);
  if (err) return 1;
//   obj_func->set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
   // creates the steepest descent optimization procedures
  SteepestDescent* pass1 = new SteepestDescent( obj_func );
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  if (err) return 1;
  pass1->set_maximum_iteration(6);
  
  QualityAssessor stop_qa=QualityAssessor(mean_ratio,QualityAssessor::MAXIMUM);
  stop_qa.add_quality_assessment(cond_num, QualityAssessor::MAXIMUM,err);
  if (err) return 1;
  
  
   //**************Set stopping criterion****************
// StoppingCriterion sc1(&stop_qa,1.0,1.8);
    //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,1);
  TerminationCriterion tc2;
  tc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err) return 1;
// CompositeAndStoppingCriterion sc(&sc1,&sc2);
  pass1->set_inner_termination_criterion(&tc2);
 // sets a culling method on the first QualityImprover
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);

  // adds 1 pass of pass1 to mesh_set1
//  queue1.add_preconditioner(pass1, err); 
//  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa,err);
  queue1.set_master_quality_improver(pass1, err); 
  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa,err);
  if (err) return 1;
  // adds 1 passes of pass2 to mesh_set1
//  mesh_set1.add_quality_pass(pass2);

  mesh->write_vtk("original_mesh", err); 
  if (err) return 1;
  
    // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err);
  if (err) return 1;
  
  mesh->write_vtk("smoothed_mesh", err); 
  if (err) return 1;
  
  return 0;
}
