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
//  LAST-MOD: 23-Jul-03 at 18:10:35 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifndef MSQ_USE_OLD_STD_HEADERS
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
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorithms
#include "ConditionNumberQualityMetric.hpp"
#include "MaxTemplate.hpp"
#include "NonSmoothSteepestDescent.hpp"

#include "MeshImpl.hpp"
using namespace Mesquite;


int main()
{     
    /* Reads a Mesh file */
  const char *file_name = 
//      "../../meshFiles/2D/VTK/equil_tri2.vtk";
//      "../../meshFiles/2D/VTK/bad_circle_tri_rhr.vtk";
//      "../../meshFiles/2D/VTK/tri_20258.vtk";
//      "../../meshFiles/3D/VTK/tet_1.vtk";
//      "../../meshFiles/3D/VTK/cube_tet_2.vtk";
     "../../meshFiles/3D/VTK/tire.vtk";
  printf("Loading mesh set 1\n");
  MsqPrintError err( cout );
  Mesquite::MeshImpl* mesh = new Mesquite::MeshImpl;
  mesh->read_vtk(file_name, err);
  if (err) return 1;
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
    //  printf("Creating mesh set 1\n");
  mesh_set1.add_mesh(mesh, err); 
  if (err) return 1;
  
    // Creates an intruction queue
    //  printf("Creating instruction queue\n");
  InstructionQueue queue1;

  // Creates a condition number quality metric 
  //  printf("Creating quality metric\n");
  ShapeQualityMetric* cond_no = new ConditionNumberQualityMetric;

  // Build an objective function with the quality metric
  //  printf("min template\n");
  MaxTemplate obj_func_min(cond_no);
  
  // Create the NonSmooth Steepest Descent procedures
  //  printf("creating optimizer\n");
  NonSmoothSteepestDescent minmax_method( &obj_func_min );
   
  minmax_method.set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err, 1);
  if (err) return 1;

  // Set a culling method on the first QualityImprover
  minmax_method.add_culling_method(PatchData::NO_BOUNDARY_VTX);

  // Set a termination criterion
  TerminationCriterion tc2;
  tc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err) return 1;
  minmax_method.set_outer_termination_criterion(&tc2);
  // Set up the quality assessor
  //  printf("Setting up the quality assessor\n");
  QualityAssessor quality_assessor=QualityAssessor(cond_no,QualityAssessor::MAXIMUM);
  quality_assessor.add_quality_assessment(cond_no,QualityAssessor::MINIMUM, err); 
  if (err) return 1;
  quality_assessor.add_quality_assessment(cond_no,QualityAssessor::AVERAGE, err);
  if (err) return 1;

  // assess the quality of the initial mesh
  queue1.add_quality_assessor(&quality_assessor, err); 
  if (err) return 1;

  // Set the max min method to be the master quality improver
  queue1.set_master_quality_improver(&minmax_method, err); 
  if (err) return 1;

  // assess the quality of the final mesh
  queue1.add_quality_assessor(&quality_assessor, err); 
  if (err) return 1;

  // write out the original mesh
  //  printf("Writing out the original mesh\n");
  mesh->write_vtk("original_mesh", err); 
  if (err) return 1;

  // launches optimization on mesh_set1
  //  printf("Running the instruction queue\n");
  queue1.run_instructions(mesh_set1, err); 
  if (err) return 1;

  // write out the smoothed mesh
  //  printf("Writing out the final mesh\n");
  mesh->write_vtk("smoothed_mesh", err); 
  if (err) return 1;
  
  delete cond_no;
  return 0;
}
