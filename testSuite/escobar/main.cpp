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
// ORIG-DATE: 11-Feb-04 at 10:57:52 by P. Knupp
//  LAST-MOD: 
//
//
// DESCRIPTION: Tests Condition Number with Escobar Variant
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifndef MSQ_USE_OLD_IO_HEADERS
#include <iostream>
using std::cout;
using std::endl;
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
#include "PlanarDomain.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "MsqError.hpp"
#include "MeshSet.hpp"
#include "ShapeImprovementWrapper.hpp"
// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"


using namespace Mesquite;

int main()
{
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  MsqPrintError err(cout);

  // If want 2D test use this section (and comment out next)

    // create geometry: plane z=0, normal (0,0,1)
    Vector3D pnt(0,0,0);
    Vector3D s_norm(0,0,1);
    Mesquite::PlanarDomain msq_geom(s_norm, pnt);
     
 //mesh->read_vtk("../../meshFiles/2D/VTK/hybrid_3quad_1tri_tangled.vtk", err);
 //mesh->read_vtk("../../meshFiles/2D/VTK/rotsq.vtk", err);
 mesh->read_vtk("../../meshFiles/2D/VTK/horseshoe.vtk", err);
 if (err) return 1;

    // initializes a MeshSet object
    MeshSet mesh_set1;
    mesh_set1.set_domain_constraint(&msq_geom, err);
    if (err) return 1;

  // End 2D Section

  // If want 3D test use this section (and comment out former)

    //mesh->read_vtk("../../meshFiles/3D/VTK/two_tets.vtk", err);
   // Note: this mesh has two inverted tets

    // initializes a MeshSet object
    //MeshSet mesh_set1;

 // End 3D Section

  mesh_set1.add_mesh(mesh, err);
  if (err) return 1;
  
    // creates an intruction queue
  InstructionQueue queue1;
  
  // Creates a condition number quality metric
  //  printf("Creating quality metric\n");
  ShapeQualityMetric* cond_no = new ConditionNumberQualityMetric;
                                                                               
    // ... and builds an objective function with it
  LPtoPTemplate obj_func(cond_no, 2, err);
  if (err) return 1;
  obj_func.set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
  
    // creates the optimization procedure
  ConjugateGradient pass1( &obj_func, err );
  //FeasibleNewton pass1( &obj_func );
  if (err) return 1;
  pass1.set_patch_type(PatchData::GLOBAL_PATCH, err);
  if (err) return 1;
  
  QualityAssessor qa=QualityAssessor(cond_no,QualityAssessor::ALL_MEASURES);
  
    // **************Set termination criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,20,err);
 
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  
  pass1.set_inner_termination_criterion(&tc_inner);
  pass1.set_outer_termination_criterion(&tc_outer);

    // sets a culling method on the first QualityImprover
    //This is an old command that still needs to be there.  It has
    //nothing to do with 'culling methods' described in TerminationCriterion.
  pass1.add_culling_method(PatchData::NO_BOUNDARY_VTX);
  queue1.add_quality_assessor(&qa,err); 
  if (err) return 1;
    // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(&pass1, err); 
  if (err) return 1;
  queue1.add_quality_assessor(&qa,err); 
  if (err) return 1;
  mesh->write_vtk("original_mesh.vtk",err); 
  if (err) return 1;
  
  queue1.run_instructions(mesh_set1, err); 
  if (err) return 1;
  mesh->write_vtk("smoothed_mesh.vtk",err); 
  if (err) return 1;

  delete cond_no;
  print_timing_diagnostics( cout );
  return 0;
}
 
