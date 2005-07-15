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
//    AUTHOR: Michael Brewer
//       ORG: Sandia National Labs
//    E-MAIL: mbrewer@sandia.gov
//
// ORIG-DATE: Jan. 29, 2003
//  LAST-MOD: 25-Feb-04 at 10:49:32 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PlanarGeometryTest.cpp

Regression testing using the planar geometry capabilities in
SimplifiedGeometryEngine.
 */
// DESCRIP-END.
//

#include "PatchDataInstances.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include <math.h>

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
//#include "StoppingCriterion.hpp"
#include "QualityAssessor.hpp"

#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "ASMQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LaplacianSmoother.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "MultiplyQualityMetric.hpp"
#include "PlanarDomain.hpp"

#include "MeshImpl.hpp"
#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else 
#include <iostream>
using std::cout;
using std::endl;
#endif

using namespace Mesquite;

class PlanarGeometryTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PlanarGeometryTest);
    //run steepest descent on the tangled_tri.vtk mesh
  CPPUNIT_TEST (test_plane_tri_tangled);
    //run cg on tangled_quad.vtk mesh
  CPPUNIT_TEST (test_plane_quad_tangled);
    //run cg with asm metric on tri mesh in y=-5 plane
  CPPUNIT_TEST (test_plane_tri_xz);
  
  CPPUNIT_TEST_SUITE_END();
  
private:
  double qualTol;//double used for double comparisons
  int pF;//PRINT_FLAG
public:
  void setUp()
  {
      //pF=1;//PRINT_FLAG IS ON
      pF=0;//PRINT_FLAG IS OFF
        //tolerance double
      qualTol=MSQ_MIN;
  }

  void tearDown()
  {
  }
  
public:
  PlanarGeometryTest()
    {}
  
   void test_plane_tri_tangled()
   {
     Mesquite::MsqPrintError err(cout); 
     Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
     
     mesh->read_vtk("../../meshFiles/2D/VTK/tangled_tri.vtk", err);
     
       // initialises a MeshSet object
     MeshSet mesh_set1;
     mesh_set1.add_mesh(mesh, err); CPPUNIT_ASSERT(!err);
     
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     
       //create geometry: plane z=5, normal (0,0,1)
     Vector3D pnt(0,0,5);
     Vector3D s_norm(0,0,1);
     Mesquite::PlanarDomain msq_geom(s_norm, pnt);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     mesh_set1.set_domain_constraint(&msq_geom, err); CPPUNIT_ASSERT(!err);
     
       // creates an intruction queue
     InstructionQueue queue1, queue2;
     
       // creates a mean ratio quality metric ...
     ShapeQualityMetric* shape = new ConditionNumberQualityMetric;
     UntangleQualityMetric* untan = new UntangleBetaQualityMetric(.1);
     
       // ... and builds an objective function with it (untangle)
     LInfTemplate* untan_func = new LInfTemplate(untan);
     untan_func->set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
     LPtoPTemplate* shape_func = new LPtoPTemplate(shape,2,err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     shape_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
       // creates the steepest descent optimization procedures
     SteepestDescent* pass1 = new SteepestDescent( untan_func );
     SteepestDescent* pass2 = new SteepestDescent( shape_func );
     pass1->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err,1 ,1);
     pass2->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     QualityAssessor stop_qa=QualityAssessor(untan,QualityAssessor::MAXIMUM, err);
     CPPUNIT_ASSERT(!err);
     QualityAssessor qa=QualityAssessor(shape,QualityAssessor::MAXIMUM, err);
     CPPUNIT_ASSERT(!err);
     if(pF==0){
       stop_qa.disable_printing_results();
       qa.disable_printing_results();
     }  
       //**********Set stopping criterion  untangle ver small ********
       //StoppingCriterion sc_qa(&stop_qa,-100,MSQ_MIN);
       //pass1->set_stopping_criterion(&sc_qa);
     TerminationCriterion sc_of;
     sc_of.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,10,err);
     pass1->set_outer_termination_criterion(&sc_of);
     
       //**********Set stopping criterion  5 iterates ****************
       //StoppingCriterion sc5(StoppingCriterion::NUMBER_OF_PASSES,5);
       //pass2->set_stopping_criterion(&sc5);
     TerminationCriterion sc5;
     sc5.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5,err);
     pass2->set_inner_termination_criterion(&sc5);
       // sets a culling method on the first QualityImprover
     pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
     pass2->add_culling_method(PatchData::NO_BOUNDARY_VTX);
       //TerminationCriterion sc_inner;
       //sc_inner.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5,err);
       //pass2->set_inner_termination_criterion(&sc_inner);
       //pass2->set_maximum_iteration(5);
  
     queue1.set_master_quality_improver(pass1, err); CPPUNIT_ASSERT(!err);
     queue2.set_master_quality_improver(pass2, err); CPPUNIT_ASSERT(!err);
       //********************UNTANGLE*******************************
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       // launches optimization on mesh_set1
     double orig_qa_val=stop_qa.loop_over_mesh(mesh_set1,err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     queue1.run_instructions(mesh_set1, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     double fin_qa_val=stop_qa.loop_over_mesh(mesh_set1,err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       //make sure 'quality' improved
     CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
       //make sure sc_qa really was satisfied
     CPPUNIT_ASSERT( fin_qa_val <= MSQ_MIN );

      //********************SMOOTH*******************************
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       // launches optimization on mesh_set1
     orig_qa_val=qa.loop_over_mesh(mesh_set1,err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     queue2.run_instructions(mesh_set1, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     fin_qa_val=qa.loop_over_mesh(mesh_set1,err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       //make sure 'quality' improved
     CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
     print_timing_diagnostics(cout);
     delete shape;
     delete untan;
     delete pass1;
     delete pass2;
     delete untan_func;
     delete shape_func;
   }
  
  void test_plane_quad_tangled()
     {
       Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
       MsqPrintError err(cout); 
       mesh->read_vtk("../../meshFiles/2D/VTK/tangled_quad.vtk", err);
       
         // initialises a MeshSet object
       MeshSet mesh_set1;
       mesh_set1.add_mesh(mesh, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //create geometry: plane z=5, normal (0,0,1)
       Vector3D pnt(0,0,5);
       Vector3D s_norm(0,0,1);
       Mesquite::PlanarDomain msq_geom(s_norm, pnt);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       mesh_set1.set_domain_constraint(&msq_geom, err); CPPUNIT_ASSERT(!err);
       
         // creates an intruction queue
       InstructionQueue queue1, queue2;

         //creates a mean ratio quality metric ...
       ShapeQualityMetric* shape = new ConditionNumberQualityMetric();
       UntangleQualityMetric* untan= new UntangleBetaQualityMetric(.1);
  
         // ... and builds an objective function with it (untangle)
       LInfTemplate* untan_func = new LInfTemplate(untan);
       untan_func->set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
       LPtoPTemplate* shape_func = new LPtoPTemplate(shape,2,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       shape_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
         // creates the cg optimization procedures
       ConjugateGradient* pass1 = new ConjugateGradient( untan_func, err );
       ConjugateGradient* pass2 = new ConjugateGradient( shape_func, err );
       pass1->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err,1 ,1);
       pass2->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       QualityAssessor stop_qa=QualityAssessor(untan,QualityAssessor::MAXIMUM, err);
       CPPUNIT_ASSERT(!err);
       QualityAssessor qa=QualityAssessor(shape,QualityAssessor::MAXIMUM, err);
       CPPUNIT_ASSERT(!err);
         //turn off printing if print flag not set.
       if(pF==0){
         stop_qa.disable_printing_results();
         qa.disable_printing_results();
       }
       //**********Set stopping criterion  untangle ver small ********
       //StoppingCriterion sc_qa(&stop_qa,-100,MSQ_MIN);
       //pass1->set_stopping_criterion(&sc_qa);
       TerminationCriterion sc_of;
       sc_of.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,10,err);
       pass1->set_outer_termination_criterion(&sc_of);
       
         //**********Set stopping criterion  5 iterates ****************
         //StoppingCriterion sc5(StoppingCriterion::NUMBER_OF_PASSES,5);
         //pass2->set_stopping_criterion(&sc5);
       TerminationCriterion sc5;
       sc5.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5,err);
       pass2->set_inner_termination_criterion(&sc5);
         // sets a culling method on the first QualityImprover
       pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
       pass2->add_culling_method(PatchData::NO_BOUNDARY_VTX);
         //pass2->set_maximum_iteration(5);
         //TerminationCriterion sc_inner;
         //sc_inner.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5,err);
         //pass2->set_inner_termination_criterion(&sc_inner);
       queue1.set_master_quality_improver(pass1, err); CPPUNIT_ASSERT(!err);
       queue2.set_master_quality_improver(pass2, err); CPPUNIT_ASSERT(!err);
         //********************UNTANGLE*******************************
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // launches optimization on mesh_set1
       double orig_qa_val=stop_qa.loop_over_mesh(mesh_set1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       queue1.run_instructions(mesh_set1, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       double fin_qa_val=stop_qa.loop_over_mesh(mesh_set1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //make sure 'quality' improved
       CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
         //make sure sc_qa really was satisfied
       CPPUNIT_ASSERT( fin_qa_val <= MSQ_MIN );
       
         //********************SMOOTH*******************************
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // launches optimization on mesh_set1
       orig_qa_val=qa.loop_over_mesh(mesh_set1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       queue2.run_instructions(mesh_set1, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
       CPPUNIT_ASSERT(!err);
       fin_qa_val=qa.loop_over_mesh(mesh_set1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //make sure 'quality' improved
       CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
       print_timing_diagnostics(cout);
       delete shape;
       delete untan;
       delete pass1;
       delete pass2;
       delete untan_func;
       delete shape_func;
     }
  
  void test_plane_tri_xz()
     {
       MsqPrintError err(cout); 
       Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
       mesh->read_vtk("../../meshFiles/2D/VTK/tri_5_xz.vtk", err);
       
         // initialises a MeshSet object
       MeshSet mesh_set1;
       mesh_set1.add_mesh(mesh, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //create geometry: plane y=5, normal (0,1,0)
       Vector3D pnt(0,-5,0);
       Vector3D s_norm(0,-1,0);
       Mesquite::PlanarDomain msq_geom(s_norm, pnt);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       mesh_set1.set_domain_constraint(&msq_geom, err); CPPUNIT_ASSERT(!err);
       
         // creates an intruction queue
       InstructionQueue queue1;
       
         //creates a asm quality metric ...
       SmoothnessQualityMetric* smooth = new ASMQualityMetric;
       
         // ... and builds an objective function with it (untangle)
       LPtoPTemplate* smooth_func = new LPtoPTemplate(smooth,1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       smooth_func->set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
         // creates the cg optimization procedures
       ConjugateGradient* pass1 = new ConjugateGradient( smooth_func, err );
         //pass1->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err,1 ,1);
       pass1->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
       pass1->set_debugging_level(1);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       QualityAssessor qa=QualityAssessor(smooth,QualityAssessor::AVERAGE, err);
       CPPUNIT_ASSERT(!err);
       
         //**********Set stopping criterion  5 iterates ****************
       TerminationCriterion sc5;
       sc5.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5,err);
       pass1->set_inner_termination_criterion(&sc5);
         //StoppingCriterion sc5(StoppingCriterion::NUMBER_OF_PASSES,5);
         //pass1->set_stopping_criterion(&sc5);
         // sets a culling method on the first QualityImprover
       pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
       TerminationCriterion sc_inner;
       sc_inner.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5,err);
       pass1->set_inner_termination_criterion(&sc_inner);
         //pass1->set_maximum_iteration(5);
       
       queue1.set_master_quality_improver(pass1, err); CPPUNIT_ASSERT(!err);
         //********************UNTANGLE*******************************
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // launches optimization on mesh_set1
       double orig_qa_val=qa.loop_over_mesh(mesh_set1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       queue1.run_instructions(mesh_set1, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       double fin_qa_val=qa.loop_over_mesh(mesh_set1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //make sure 'quality' improved
       CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
       print_timing_diagnostics(cout);
       delete smooth;
       delete smooth_func;
       delete pass1;
     }
  
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PlanarGeometryTest, "PlanarGeometryTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PlanarGeometryTest, "Regression");
