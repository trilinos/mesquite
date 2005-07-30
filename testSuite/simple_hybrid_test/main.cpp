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
//  LAST-MOD: 23-Jul-03 at 18:11:05 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "PlanarDomain.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "MsqError.hpp"
#include "ShapeImprovementWrapper.hpp"
// algorythms
#include "IdealWeightInverseMeanRatio.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"


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

using namespace Mesquite;

int main()
{
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  MsqPrintError err(cout);
    //create geometry: plane z=0, normal (0,0,1)
  Vector3D pnt(0,0,0);
  Vector3D s_norm(0,0,1);
  Mesquite::PlanarDomain msq_geom(s_norm, pnt);
     
    //mesh->read_vtk("../../meshFiles/2D/VTK/cube-clip-corner.vtk", err);
  mesh->read_vtk("../../meshFiles/2D/VTK/hybrid_3quad_1tri.vtk", err);
  if (err) return 1;
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
    //   mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    //   mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
    //mean_ratio->set_averaging_method(QualityMetric::SUM, err);
    //MSQ_CHKERR(err);
  
    // ... and builds an objective function with it
  LPtoPTemplate obj_func(mean_ratio, 2, err);
  if (err) return 1;;
  obj_func.set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  
    // creates the steepest descent, feas newt optimization procedures
    //ConjugateGradient* pass1 = new ConjugateGradient( &obj_func, err );
  FeasibleNewton pass1( &obj_func );
  pass1.set_patch_type(PatchData::GLOBAL_PATCH, err);
  if (err) return 1;;
  
  QualityAssessor qa=QualityAssessor(mean_ratio,QualityAssessor::ALL_MEASURES, err);
  if (err) return 1;;
  
    // **************Set termination criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err) return 1;;
    //_inner.add_criterion_type_with_double(
    // TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, OF_value, err);
    //tc_inner.add_criterion_type_with_double(
    //TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE, OF_value, err);
  TerminationCriterion tc_outer;
    //tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err) return 1;;
  
  pass1.set_inner_termination_criterion(&tc_inner);
  pass1.set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&qa,err); 
  if (err) return 1;
    // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(&pass1, err);
  if (err) return 1;
  queue1.add_quality_assessor(&qa,err); 
  if (err) return 1;
  mesh->write_vtk("original_mesh.vtk",err); 
  if (err) return 1;
  
  queue1.run_instructions(mesh, &msq_geom, err); 
  if (err) return 1;
  mesh->write_vtk("smoothed_mesh.vtk",err); 
  if (err) return 1;
    //std::cout<<"\n\nNow running the shape wrapper.\n=n";
    //ShapeImprovementWrapper wrap(100);
    //wrap.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  delete mean_ratio;
  print_timing_diagnostics(cout);
  return 0;
}
 
