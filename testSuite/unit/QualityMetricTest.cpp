// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Michael Brewer
//       ORG: Sandia National Labs
//    E-MAIL: mbrewer@sandia.gov
//
// ORIG-DATE: 03-Dec-02
//  LAST-MOD: 03-Dec-02 
//
// DESCRIPTION:
// ============
/*! \file QualityMetricTest.cpp

Unit testing of various QualityMetrics primarily to test for
correct metric return values. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "QualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "GeneralizedConditionNumberQualityMetric.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "InverseMeanRatioQualityMetric.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"
#include <math.h>

using namespace Mesquite;

class QualityMetricTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(QualityMetricTest);
  CPPUNIT_TEST (test_condition_number);
  CPPUNIT_TEST_SUITE_END();
  
private:
  
  PatchData triPatch;
  PatchData quadPatch;
  
  
   


public:
  void setUp()
  {
    MsqError err;
    int ind[20];
    
    int elem_ind[8];
    
     /* Our triangular patch is made of two tris.  triE1 is a perfect
        equilateral (the ideal for most metrics).  triE2 is an arbitrary
        triangle for which we calculate the metric value.
     */
     triPatch.reserve_vertex_capacity(4, err); MSQ_CHKERR(err);
     ind[0]=triPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=triPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=triPatch.add_vertex(NULL, NULL, .5,sqrt(3.0)/2.0,0.0,
                                true, err);
     MSQ_CHKERR(err);
     ind[3]=triPatch.add_vertex(NULL, NULL, 2.0,2.0,2.0, true, err);
     MSQ_CHKERR(err);
     triPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal equilateral
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     triPatch.add_element(NULL, NULL, elem_ind, TRIANGLE, err);
     MSQ_CHKERR(err);
       //add "arbitrary" tri
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[3];
     elem_ind[2] = ind[1];
     triPatch.add_element(NULL, NULL, elem_ind, TRIANGLE, err);
     MSQ_CHKERR(err);
  }

  void tearDown()
  {
  }
  
public:
  QualityMetricTest()
    {}
  
   void test_condition_number()
   {
      MsqError err;
      double val;
      MsqMeshEntity* elems;
      MsqVertex* verts = triPatch.get_vertex_array(err);
      elems=triPatch.get_element_array(err);
      MSQ_CHKERR(err);
      ShapeQualityMetric *met = ConditionNumberQualityMetric::create_new();
      ShapeQualityMetric *gmet = GeneralizedConditionNumberQualityMetric::create_new();
        //Check condition number of ideal tri
      val = met->evaluate_element(triPatch,&elems[0],err);MSQ_CHKERR(err);
      CPPUNIT_ASSERT(fabs(val-1.0)<MSQ_MIN);
        //For now, make sure cond num and generalized cond num give
        //equivalent answer for arbitrary tri.
      val = met->evaluate_element(triPatch,&elems[1],err); MSQ_CHKERR(err);
      val -= gmet->evaluate_element(triPatch,&elems[1],err); MSQ_CHKERR(err);
      CPPUNIT_ASSERT(fabs(val)<MSQ_MIN);
   }
  
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "QualityMetricTest");

