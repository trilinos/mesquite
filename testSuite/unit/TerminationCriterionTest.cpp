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
//  LAST-MOD: March 5, 2003
//
// DESCRIPTION:
// ============
/*! \file TerminationCriterionTest.cpp

Tests for the TerminationCriterion class.. 

*/
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
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "ConjugateGradient.hpp"
#include "SimplifiedGeometryEngine.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"
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
      char file_name[128];
        /*Reads a TSTT Mesh file */
      TSTT::Mesh_Handle mesh;
      TSTT::MeshError tstt_err;
      TSTT::Mesh_Create(&mesh, &tstt_err);
      strcpy(file_name, "../../meshFiles/2D/VTK/tri_5_xz.vtk");
      TSTT::Mesh_Load(mesh, file_name, &tstt_err);
  
        // initialises a MeshSet object
      MeshSet mesh_set1;
      mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
        //create geometry
      Vector3D center(2,2,0);
      Vector3D geo_center(0,0,0);
      Vector3D pnt(0,-5,0);
      Vector3D s_norm(0, -1,0);
      SimplifiedGeometryEngine msq_geom;
      msq_geom.set_geometry_to_plane(s_norm,pnt,err);
      mesh_set1.set_simplified_geometry_engine(&msq_geom);
  
        // creates an intruction queue        
      InstructionQueue queue1;
      
        //     creates a mean ratio quality metric ...
      ShapeQualityMetric* cond_num=ConditionNumberQualityMetric::create_new();
      LPTemplate* obj_func = new LPTemplate(cond_num, 2, err);
      CPPUNIT_ASSERT(!err.errorOn);
      obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
      ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
      CPPUNIT_ASSERT(!err.errorOn);
      pass1->set_outer_termination_criterion(tc_outer);
      pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
      
      pass1->set_debugging_level(0);
      QualityAssessor all_qa=QualityAssessor(cond_num,
                                             QualityAssessor::ALL_MEASURES);
      queue1.add_quality_assessor(&all_qa,err);
      CPPUNIT_ASSERT(!err.errorOn);
      queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
      CPPUNIT_ASSERT(!err.errorOn);
        //queue1.add_quality_assessor(&all_qa,err);
      CPPUNIT_ASSERT(!err.errorOn);
      queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
      CPPUNIT_ASSERT(!err.errorOn);
      delete pass1;
      delete obj_func;
      delete cond_num;
    }
  
    //NUMBER OF ITERATES
  void test_number_of_iterates()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,
                                     2, err);
      if(pF)
        std::cout<<"\nTEST_NUMBER_OF_ITERATES\n";
      test_outer_criterion(&t1,err);
    }
  
    //GRADIENT NORM ABSOLUTE
  void test_gradient_norm_absolute()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::GRADIENT_NORM_ABSOLUTE,.02, err);
      if(pF)
        std::cout<<"\nTEST_GRADIENT_NORM_ABS)LUTE\n";
      test_outer_criterion(&t1,err);
    }
  
    //GRADIENT NORM RELATIVE
  void test_gradient_norm_relative()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::GRADIENT_NORM_RELATIVE,.3, err);
      if(pF)
        std::cout<<"\nTEST_GRADIENT_NORM_RELATIVE\n";
      test_outer_criterion(&t1,err);
    }
  
    //CPU TIME
  void test_cpu_time()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::CPU_TIME,5, err);
      if(pF)
        std::cout<<"\nTEST_CPU_TIME\n";
      test_outer_criterion(&t1,err);
    }
  
    //CPU TIME OR NUMBER OF ITERATES
   void test_cpu_time_or_iterates()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::CPU_TIME,5, err);
      t1.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5000, err);
      if(pF)
        std::cout<<"\nTEST_CPU_TIME_OR_ITERATES\n";
      test_outer_criterion(&t1,err);
    }
  
    //NUMBER OF ITERATES OR CPU TIME
  void test_iterates_or_cpu_time()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::CPU_TIME,5000, err);
      t1.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,2, err);
      if(pF)
        std::cout<<"\nTEST_ITERATES_OR_CPU_TIME\n";
      test_outer_criterion(&t1,err);
    }
  
    //QUALITY IMPROVEMENT ABSOLUTE
  void test_quality_improvement_absolute()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE,11.95, err);
      if(pF)
        std::cout<<"\nTEST_QUALITY_IMPROVEMENT_ABSOLUTE\n";
      test_outer_criterion(&t1,err);
    }

    //QUALITY IMPROVEMENT RELATIVE
  void test_quality_improvement_relative()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_RELATIVE,.996, err);
      if(pF)
        std::cout<<"\nTEST_QUALITY_IMPROVEMENT_RELATIVE\n";
      test_outer_criterion(&t1,err);
    }

    //SUCCESSIVE IMPROVEMENTS ABSOLUTE
  void test_successive_improvements_absolute()
    {
      MsqError err;
      TerminationCriterion t1;
      t1.add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_ABSOLUTE,.005, err);
      if(pF)
        std::cout<<"\nTEST_SUCCESSIVE_IMPROVEMENTS_ABSOLUTE\n";
      test_outer_criterion(&t1,err);
    }
  
    //SUCCESSIVE IMPROVEMENTS RELATIVE
  void test_successive_improvements_relative()
     {
       MsqError err;
       TerminationCriterion t1;
       t1.add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE,.25, err);
       if(pF)
         std::cout<<"\nTEST_SUCCESSIVE_IMPROVEMENTS_RELATIVE\n";
       test_outer_criterion(&t1,err);
     }
  
    //VERTEX BOUND OR CPU TIME
  void test_vertex_bound_or_cpu_time()
     {
       MsqError err;
       TerminationCriterion t1;
       t1.add_criterion_type_with_double(TerminationCriterion::CPU_TIME,5000, err);
       t1.add_criterion_type_with_double(TerminationCriterion::BOUNDED_VERTEX_MOVEMENT, .01 , err);
       if(pF)
        std::cout<<"\nTEST_VERTEX_BOUND_OR_CPU_TIME\n";
       test_outer_criterion(&t1,err);
     }
    //CPU TIME OR VERTEX BOUND
  void test_cpu_time_or_vertex_bound()
     {
       MsqError err;
       TerminationCriterion t1;
       t1.add_criterion_type_with_double(TerminationCriterion::CPU_TIME,5, err);
       t1.add_criterion_type_with_double(TerminationCriterion::BOUNDED_VERTEX_MOVEMENT, 50000 , err);
       if(pF)
         std::cout<<"\nTEST_CPU_TIME_OR_VERTEX_BOUND\n";
       test_outer_criterion(&t1,err);
     }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TerminationCriterionTest, "TerminationCriterionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TerminationCriterionTest, "Regression");
