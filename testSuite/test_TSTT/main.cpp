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
//  LAST-MOD: 10-Feb-04 at 22:44:58 by Thomas Leurent
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


#ifdef MSQ_USE_TSTT_OVERTURE_IMPL
#include "TSTT_Overture_Mesh.hh"
#endif

#ifdef MSQ_USE_TSTT_AOMD_IMPL
#include "TSTT_LocalTSTTMesh_Impl.hh"
#endif

#include "Mesquite.hpp"
#include "MeshTSTT.hpp"
#include "MesquiteError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"

// algorythms
#include "MeanRatioQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"

#include "MsqMessage.hpp"


using namespace Mesquite;

using std::cout;
using std::endl;

#undef __FUNC__
#define __FUNC__ "main"
int main(int argc, char* argv[])
{
  Mesquite::MsqError err;
  char file_name[256];
  double OF_value = 1.;
  
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
  
  TSTT::Mesh tstt_mesh;

#ifdef MSQ_USE_TSTT_OVERTURE_IMPL
  TSTT_Overture::Mesh ov_mesh = TSTT_Overture::Mesh::_create();
  tstt_mesh = ov_mesh;
#endif

#ifdef MSQ_USE_TSTT_AOMD_IMPL
  TSTT::LocalTSTTMesh aomd_mesh = TSTT::LocalTSTTMesh::_create();
  tstt_mesh = aomd_mesh;
#endif
  
  TSTT::CoreEntitySetQuery tstt_core_query = tstt_mesh;

  if ( !tstt_core_query )
  {
     err.set_msg("ERROR: could not perform tstt cast to CoreEntitySetQuery");
     return 1;
  }

  tstt_core_query.load(file_name);

  Mesquite::MeshTSTT* mesh = new Mesquite::MeshTSTT(tstt_mesh, err);
  if (err.errorOn) return 1;;
  
  // Creates a domain constraint
  Vector3D normal(0,0,1);
  Vector3D point(0,0,0);
  Mesquite::PlanarDomain planar_domain(normal,point,mesh);

  // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); 
  if (err.errorOn) return 1;;
  //  mesh_set1.set_domain_constraint(&planar_domain, err); MSQ_CHKERR(err);

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
//   SmoothnessQualityMetric* mean_ratio = new EdgeLengthQualityMetric;
  ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
//  mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
//   mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
  mean_ratio->set_averaging_method(QualityMetric::SUM, err); 
  if (err.errorOn) return 1;;
  
  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err.errorOn) return 1;;
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  
  // creates the steepest descentfeas newt optimization procedures
//  ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  if (err.errorOn) return 1;;
  
  QualityAssessor stop_qa=QualityAssessor(mean_ratio,QualityAssessor::AVERAGE);
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_double(
           TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, OF_value, err);
  if (err.errorOn) return 1;;
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err.errorOn) return 1;;
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  // sets a culling method on the first QualityImprover
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  
  queue1.add_quality_assessor(&stop_qa, err);
  if (err.errorOn) return 1;;
   
  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err);
  if (err.errorOn) return 1;;
  
  queue1.add_quality_assessor(&stop_qa, err); 
  if (err.errorOn) return 1;;

  mesh_set1.write_vtk("original", err); 
  if (err.errorOn) return 1;;

  // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err);
  if (err.errorOn) return 1;;

  mesh_set1.write_vtk("smoothed", err); 
  if (err.errorOn) return 1;;

  Message::print_timing_diagnostics();
  return 0;
}
