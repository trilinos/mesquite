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
//    AUTHOR: Michael Brewer
//       ORG: Sanida National Laboratories
//    E-MAIL: mbrewer@sandia.gov
//
// ORIG-DATE: March 5, 2003
//  LAST-MOD: 25-Feb-04 at 10:50:02 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file TerminationCriterionTest.cpp

Tests for the TerminationCriterion class.. 

*/

#include "meshfiles.h"

//
#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "MeshImpl.hpp"

// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "ConjugateGradient.hpp"
#include "PlanarDomain.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "MsqFreeVertexIndexIterator.hpp"
#include <list>
#include <iterator>

using namespace Mesquite;
class TerminationCriterionTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(TerminationCriterionTest);
    //simple test with number of iterates = 2 termination criterion
  CPPUNIT_TEST (test_number_of_iterates);
    //terminate when the norm of the gradient is below .02
  CPPUNIT_TEST (test_gradient_norm_absolute);
    //terminate when the norm of the gradient is below .3*initialNorm
  CPPUNIT_TEST (test_gradient_norm_relative);
    //terminate when the time exceeds 5 seconds
  CPPUNIT_TEST (test_cpu_time);
    //terminate when the time exceeds 5 seconds or iterates exceed 5000
  CPPUNIT_TEST (test_cpu_time_or_iterates);
    //terminate when the iterates exceed 2 or the time exceeds 5000 seconds
  CPPUNIT_TEST (test_iterates_or_cpu_time);
    //terminate when the function value is below 11.95
  CPPUNIT_TEST (test_quality_improvement_absolute);
    //terminate when the function value is below (.995*intialOFValue)
  CPPUNIT_TEST (test_quality_improvement_relative);
    //terminate when the successive improvements is below .005
  CPPUNIT_TEST (test_successive_improvements_absolute);
    //terminat when the successive imporovements is below .996 of the original
  CPPUNIT_TEST (test_successive_improvements_relative);
    //tests vertex bound should NOT stop due to vertices outside the bound
  CPPUNIT_TEST (test_cpu_time_or_vertex_bound);
    //tests vertex bound SHOULD stop due to vertices outside the bound
  CPPUNIT_TEST (test_vertex_bound_or_cpu_time);
  CPPUNIT_TEST_SUITE_END();
  int pF;
  
public:
  void setUp()
  {
    pF=1;//PRINT_FLAG IS ON
      //pF=0;//PRINT_FLAG IS OFF
  }

  void tearDown()
  {
  }
    //Function that gets passed a pointer to the outer termination criterion.
    //This function is called by several of the tests in this suite.
  void test_outer_criterion(TerminationCriterion* tc_outer, MsqError &err)
    {
      Mesquite::MeshImpl mesh;
      mesh.read_vtk(MESH_FILES_DIR "2D/VTK/tri_5_xz.vtk", err);
      
      Vector3D pnt(0,-5,0);
      Vector3D s_norm(0, -1,0);
      Mesquite::PlanarDomain msq_geom(s_norm, pnt);
      
        // create an intruction queue        
      InstructionQueue queue1;
      
        // create a mean ratio quality metric ...
      ConditionNumberQualityMetric cond_num;
      LPtoPTemplate obj_func(&cond_num, 2, err);
      CPPUNIT_ASSERT(!err);
      ConjugateGradient pass1( &obj_func, err );
      CPPUNIT_ASSERT(!err);
      pass1.set_outer_termination_criterion(tc_outer);
      
      pass1.set_debugging_level(0);
      QualityAssessor all_qa=QualityAssessor( &cond_num );
      queue1.add_quality_assessor(&all_qa,err); CPPUNIT_ASSERT(!err);
      queue1.set_master_quality_improver(&pass1, err); CPPUNIT_ASSERT(!err);
        //queue1.add_quality_assessor(&all_qa,err); CPPUNIT_ASSERT(!err);
      queue1.run_instructions(&mesh, &msq_geom, err); CPPUNIT_ASSERT(!err);
    }
  
    //NUMBER OF ITERATES
  void test_number_of_iterates()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_iteration_limit( 2 );
      if(pF)
        std::cout<<"\nTEST_NUMBER_OF_ITERATES\n";
      test_outer_criterion(&t1,err);CPPUNIT_ASSERT(!err);
    }
  
    //GRADIENT NORM ABSOLUTE
  void test_gradient_norm_absolute()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_absolute_gradient_inf_norm( 0.02 );
      if(pF)
        std::cout<<"\nTEST_GRADIENT_INF_NORM_ABSOLUTE\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
    }
  
    //GRADIENT NORM RELATIVE
  void test_gradient_norm_relative()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_relative_gradient_inf_norm( 0.3 );
      if(pF)
        std::cout<<"\nTEST_GRADIENT_INF_NORM_RELATIVE\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
    }
  
    //CPU TIME
  void test_cpu_time()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_cpu_time( 5 );
      if(pF)
        std::cout<<"\nTEST_CPU_TIME\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
    }
  
    //CPU TIME OR NUMBER OF ITERATES
   void test_cpu_time_or_iterates()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_cpu_time( 5 );
      t1.add_iteration_limit( 5000 );
      if(pF)
        std::cout<<"\nTEST_CPU_TIME_OR_ITERATES\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
    }
  
    //NUMBER OF ITERATES OR CPU TIME
  void test_iterates_or_cpu_time()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_cpu_time( 5000 );
      t1.add_iteration_limit( 2 );
      if(pF)
        std::cout<<"\nTEST_ITERATES_OR_CPU_TIME\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
    }
  
    //QUALITY IMPROVEMENT ABSOLUTE
  void test_quality_improvement_absolute()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_absolute_quality_improvement( 143.0 );
      if(pF)
        std::cout<<"\nTEST_QUALITY_IMPROVEMENT_ABSOLUTE\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
    }

    //QUALITY IMPROVEMENT RELATIVE
  void test_quality_improvement_relative()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_relative_quality_improvement( 0.996 );
      if(pF)
        std::cout<<"\nTEST_QUALITY_IMPROVEMENT_RELATIVE\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
    }

    //SUCCESSIVE IMPROVEMENTS ABSOLUTE
  void test_successive_improvements_absolute()
    {
      MsqPrintError err(cout);
      TerminationCriterion t1;
      t1.add_absolute_successive_improvement( 0.005 );
      if(pF)
        std::cout<<"\nTEST_SUCCESSIVE_IMPROVEMENTS_ABSOLUTE\n";
      test_outer_criterion(&t1,err);
      CPPUNIT_ASSERT(!err);
   }
  
    //SUCCESSIVE IMPROVEMENTS RELATIVE
  void test_successive_improvements_relative()
     {
       MsqPrintError err(cout);
       TerminationCriterion t1;
       t1.add_relative_successive_improvement( 0.25 );
       if(pF)
         std::cout<<"\nTEST_SUCCESSIVE_IMPROVEMENTS_RELATIVE\n";
       test_outer_criterion(&t1,err);
       CPPUNIT_ASSERT(!err);
     }
  
    //VERTEX BOUND OR CPU TIME
  void test_vertex_bound_or_cpu_time()
     {
       MsqPrintError err(cout);
       TerminationCriterion t1;
       t1.add_cpu_time( 5000 );
       t1.add_bounded_vertex_movement( 0.01 );
       if(pF)
        std::cout<<"\nTEST_VERTEX_BOUND_OR_CPU_TIME\n";
       test_outer_criterion(&t1,err);
       CPPUNIT_ASSERT(!err);
     }
    //CPU TIME OR VERTEX BOUND
  void test_cpu_time_or_vertex_bound()
     {
       MsqPrintError err(cout);
       TerminationCriterion t1;
       t1.add_cpu_time( 5 );
       t1.add_bounded_vertex_movement( 50000 );
       if(pF)
         std::cout<<"\nTEST_CPU_TIME_OR_VERTEX_BOUND\n";
       test_outer_criterion(&t1,err);
       CPPUNIT_ASSERT(!err);
     }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TerminationCriterionTest, "TerminationCriterionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TerminationCriterionTest, "Regression");
