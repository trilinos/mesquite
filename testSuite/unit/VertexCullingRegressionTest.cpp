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
// ORIG-DATE: May 8, 2003
//  LAST-MOD: 23-Jul-03 at 17:44:57 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file VertexCullingRegressionTest.cpp

Regression testing using the vertex culling algorithms. 
 */
// DESCRIP-END.
//

#include "PatchDataInstances.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include <math.h>

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "QualityAssessor.hpp"

#include "EdgeLengthQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "LaplacianSmoother.hpp"
#include "PlanarDomain.hpp"
#include "TerminationCriterion.hpp"
#include "MeshImpl.hpp"

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else 
#include <iostream>
using std::cout;
using std::endl;
#endif


using namespace Mesquite;

class VertexCullingRegressionTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(VertexCullingRegressionTest);
  CPPUNIT_TEST (test_laplacian_smoothing_with_cull);
  CPPUNIT_TEST_SUITE_END();
  
private:
  int pF;//PRINT_FLAG
public:
  void setUp()
  {
      //pF=1;//PRINT_FLAG IS ON
      pF=0;//PRINT_FLAG IS OFF
  }

  void tearDown()
  {
  }
  
public:
  VertexCullingRegressionTest()
    {}

  void test_laplacian_smoothing_with_cull()
    {
        /* Read a VTK Mesh file */
      MsqPrintError err(cout);
      Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
      mesh->read_vtk("../../meshFiles/2D/VTK/square_quad_10_rand.vtk", err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err);
        // initialises a MeshSet object
      MeshSet mesh_set1;
      mesh_set1.add_mesh(mesh, err); CPPUNIT_ASSERT(!err);
      
        // creates an intruction queue
      InstructionQueue queue1;
      
        // creates a mean ratio quality metric ...
      ShapeQualityMetric* shape_metric = new ConditionNumberQualityMetric;
      SmoothnessQualityMetric* lapl_met = new EdgeLengthQualityMetric;
      lapl_met->set_averaging_method(QualityMetric::RMS,err); CPPUNIT_ASSERT(!err);
      
        // creates the laplacian smoother  procedures
      LaplacianSmoother* lapl1 = new LaplacianSmoother(err); CPPUNIT_ASSERT(!err);
      LaplacianSmoother* lapl2 = new LaplacianSmoother(err); CPPUNIT_ASSERT(!err);
      QualityAssessor stop_qa=QualityAssessor(shape_metric,QualityAssessor::MAXIMUM); CPPUNIT_ASSERT(!err);
      stop_qa.add_quality_assessment(lapl_met,QualityAssessor::ALL_MEASURES,err); CPPUNIT_ASSERT(!err);
      
        //**************Set termination criterion****************
      TerminationCriterion sc2;
      sc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1000,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err);
      sc2.add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_ABSOLUTE,0.0,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err);
        //set a criterion with a culling method for the inner criterion
      TerminationCriterion sc_cull;
      sc_cull.set_culling_type(TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE,.1,err);
      CPPUNIT_ASSERT(!err);
      TerminationCriterion sc_cull_2;
      sc_cull_2.set_culling_type(TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE,.000001,err);
      CPPUNIT_ASSERT(!err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err);
      lapl1->set_outer_termination_criterion(&sc2);
      lapl2->set_outer_termination_criterion(&sc2);
      lapl1->set_inner_termination_criterion(&sc_cull);
      lapl2->set_inner_termination_criterion(&sc_cull_2);
        // sets a culling method on the first QualityImprover
      lapl1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
      lapl2->add_culling_method(PatchData::NO_BOUNDARY_VTX);
        // adds 1 pass of pass1 to mesh_set1
      queue1.add_quality_assessor(&stop_qa,err);
       //Make sure no errors
      CPPUNIT_ASSERT(!err);
      queue1.add_preconditioner(lapl1,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err);
      queue1.set_master_quality_improver(lapl2, err);
       //Make sure no errors
      CPPUNIT_ASSERT(!err);
      queue1.add_quality_assessor(&stop_qa,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err);
      queue1.run_instructions(mesh_set1, err);
      CPPUNIT_ASSERT(!err);
      
        //Make sure no errors
      CPPUNIT_ASSERT(!err);
      delete shape_metric;
      delete lapl1;
      delete lapl2;
      delete lapl_met;
    }
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VertexCullingRegressionTest, "VertexCullingRegressionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VertexCullingRegressionTest, "Regression");
