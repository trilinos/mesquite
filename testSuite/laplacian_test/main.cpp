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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 23-Jul-03 at 18:04:37 by Thomas Leurent
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
using std::cout;
using std::endl;
#else
#include <iostream.h>
#endif

#ifdef MSQ_USE_OLD_C_HEADERS
#include <cstdlib>
#else
#include <stdlib.h>
#endif


#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "MeshImpl.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "LaplacianSmoother.hpp"
#include "EdgeLengthQualityMetric.hpp"
using namespace Mesquite;


int main()
{     
    /* Read a VTK Mesh file */
  MsqPrintError err(cout);
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk("../../meshFiles/2D/VTK/square_quad_2.vtk", err);
  if (err) return 1;
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); 
  if (err) return 1;
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ShapeQualityMetric* shape_metric = new ConditionNumberQualityMetric;
  SmoothnessQualityMetric* lapl_met = new EdgeLengthQualityMetric;
  lapl_met->set_averaging_method(QualityMetric::RMS,err);
   if (err) return 1;
 
    // creates the laplacian smoother  procedures
  LaplacianSmoother lapl1(err);
  if (err) return 1;
  QualityAssessor stop_qa=QualityAssessor(shape_metric,QualityAssessor::MAXIMUM, err);
  if (err) return 1;
  stop_qa.add_quality_assessment(lapl_met,QualityAssessor::ALL_MEASURES,err);
  if (err) return 1;
  
    //**************Set stopping criterion****************
  TerminationCriterion sc2;
  sc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,10,err);
  if (err) return 1;
  lapl1.set_outer_termination_criterion(&sc2);
    // sets a culling method on the first QualityImprover
  lapl1.add_culling_method(PatchData::NO_BOUNDARY_VTX);
  
    // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
  queue1.set_master_quality_improver(&lapl1, err); 
  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
    // adds 1 passes of pass2 to mesh_set1
    //  mesh_set1.add_quality_pass(pass2);
  
    //writeVtkMesh("original_mesh", mesh, err); MSQ_CHKERR(err);
  
    // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); 
  if (err) return 1;
  
  mesh->write_vtk("smoothed_mesh.vtk", err); 
  if (err) return 1;
  
  return 0;
}
