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
    //Test condition number and generalized condition number metrics
  CPPUNIT_TEST (test_condition_number);
    //Test mean ratio and inverse mean ratio metrics
  CPPUNIT_TEST (test_mean_ratio);
    //Test apect ratio gamma (Tri's and Tet's only)
  CPPUNIT_TEST (test_aspect_ratio_gamma);
  
  CPPUNIT_TEST_SUITE_END();
  
private:
  
  PatchData triPatch;
  PatchData quadPatch;
  PatchData tetPatch;
  PatchData hexPatch;
    //Tol used for double comparisons
  double qualTol;
public:
  void setUp()
  {
    
    qualTol = MSQ_MIN;
    MsqError err;
    int ind[20];
    
    int elem_ind[8];
    
     /* Our triangular patch is made of two tris.  tri_1 is a perfect
        equilateral (the ideal for most metrics).  tri_2 is an arbitrary
        triangle.
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

     /* Our quad patch is made of two quads.  quad_1 is a perfect
        square (the ideal for most metrics).  quad_2 is an arbitrary
        quad.
     */
     quadPatch.reserve_vertex_capacity(6, err); MSQ_CHKERR(err);
     ind[0]=quadPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=quadPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=quadPatch.add_vertex(NULL, NULL, 1.0,1.0,0.0,true, err);
     MSQ_CHKERR(err);
     ind[3]=quadPatch.add_vertex(NULL, NULL, 0.0,1.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[4]=quadPatch.add_vertex(NULL, NULL, 2.0,-1.0,.5, true, err);
     MSQ_CHKERR(err);
     ind[5]=quadPatch.add_vertex(NULL, NULL, 1.5,1.0,1.0, true, err);
     MSQ_CHKERR(err);
     quadPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal quad
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     
     quadPatch.add_element(NULL, NULL, elem_ind, QUADRILATERAL, err);
     MSQ_CHKERR(err);
       //add "arbitrary" quad
     elem_ind[0] = ind[1];
     elem_ind[1] = ind[4];
     elem_ind[2] = ind[5];
     elem_ind[3] = ind[2];
     
     quadPatch.add_element(NULL, NULL, elem_ind, QUADRILATERAL, err);
     MSQ_CHKERR(err);

     /* Our tet patch is made of two tets.  tet_1 is a perfect
        equilateral (the ideal for most metrics).  tet_2 is an arbitrary
        tet.
     */
     tetPatch.reserve_vertex_capacity(5, err); MSQ_CHKERR(err);
     ind[0]=tetPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=tetPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=tetPatch.add_vertex(NULL, NULL, 0.5,sqrt(3.0)/2.0,0.0,true, err);
     MSQ_CHKERR(err);
     ind[3]=tetPatch.add_vertex(NULL, NULL, .5, sqrt(3.0)/6.0,
                                 sqrt(2.0)/sqrt(3.0), true, err);
     MSQ_CHKERR(err);
     ind[4]=tetPatch.add_vertex(NULL, NULL, 2.0,3.0,-.5, true, err);
     MSQ_CHKERR(err);
    
     tetPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal tet
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     
     tetPatch.add_element(NULL, NULL, elem_ind, TETRAHEDRON, err);
     MSQ_CHKERR(err);
       //add "arbitrary" tet
     elem_ind[0] = ind[1];
     elem_ind[1] = ind[4];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     
     tetPatch.add_element(NULL, NULL, elem_ind, TETRAHEDRON, err);
     MSQ_CHKERR(err);

     /* Our hex patch is made of two hexes.  hex_1 is a perfect
        unit cube (the ideal for most metrics).  hex_2 is an arbitrary
        tet.
     */
     hexPatch.reserve_vertex_capacity(12, err); MSQ_CHKERR(err);
     ind[0]=hexPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=hexPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=hexPatch.add_vertex(NULL, NULL, 1.0, 1.0, 0.0, true, err);
     MSQ_CHKERR(err);
     ind[3]=hexPatch.add_vertex(NULL, NULL, 0.0, 1.0, 0.0, true, err);
     MSQ_CHKERR(err);
     ind[4]=hexPatch.add_vertex(NULL, NULL, 0.0,0.0,1.0, true, err);
     MSQ_CHKERR(err);
     ind[5]=hexPatch.add_vertex(NULL, NULL, 1.0,0.0,1.0, true, err);
     MSQ_CHKERR(err);
     ind[6]=hexPatch.add_vertex(NULL, NULL, 1.0, 1.0, 1.0, true, err);
     MSQ_CHKERR(err);
     ind[7]=hexPatch.add_vertex(NULL, NULL, 0.0, 1.0, 1.0, true, err);
     MSQ_CHKERR(err);
     ind[8]=hexPatch.add_vertex(NULL, NULL, 2.0,0.0,-1.0, true, err);
     MSQ_CHKERR(err);
     ind[9]=hexPatch.add_vertex(NULL, NULL, 3.0,2.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[10]=hexPatch.add_vertex(NULL, NULL, 2.0, 2.0, 2.0, true, err);
     MSQ_CHKERR(err);
     ind[11]=hexPatch.add_vertex(NULL, NULL, 3.0, 0.0, 3.0, true, err);
     MSQ_CHKERR(err);
    
     hexPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal hex
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     elem_ind[4] = ind[4];
     elem_ind[5] = ind[5];
     elem_ind[6] = ind[6];
     elem_ind[7] = ind[7];
     
     hexPatch.add_element(NULL, NULL, elem_ind, HEXAHEDRON, err);
     MSQ_CHKERR(err);
       //add "arbitrary" hex
     elem_ind[0] = ind[1];
     elem_ind[1] = ind[8];
     elem_ind[2] = ind[9];
     elem_ind[3] = ind[2];
     elem_ind[4] = ind[5];
     elem_ind[5] = ind[10];
     elem_ind[6] = ind[11];
     elem_ind[7] = ind[6];
     
     hexPatch.add_element(NULL, NULL, elem_ind, HEXAHEDRON, err);
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
       //START WITH TRI's
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
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tri.
     val = met->evaluate_element(triPatch,&elems[1],err); MSQ_CHKERR(err);
     val -= gmet->evaluate_element(triPatch,&elems[1],err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val)<qualTol);
     
       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check condition number of ideal quad
     val = met->evaluate_element(quadPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary quad.
     val = met->evaluate_element(quadPatch,&elems[1],err); MSQ_CHKERR(err);
     val -= gmet->evaluate_element(quadPatch,&elems[1],err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val)<qualTol);

       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check condition number of ideal tet
     val = met->evaluate_element(tetPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tet.
     val = met->evaluate_element(tetPatch,&elems[1],err); MSQ_CHKERR(err);
     val -= gmet->evaluate_element(tetPatch,&elems[1],err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val)<qualTol);

       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check condition number of ideal hex
     val = met->evaluate_element(hexPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tet.
     val = met->evaluate_element(hexPatch,&elems[1],err); MSQ_CHKERR(err);
     val -= gmet->evaluate_element(hexPatch,&elems[1],err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val)<qualTol);
   }


   void test_mean_ratio()
   {
       //START WITH TRI's
     MsqError err;
     double val;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = MeanRatioQualityMetric::create_new();
     ShapeQualityMetric *imet = InverseMeanRatioQualityMetric::create_new();
       //Check mean ratio of ideal tri
     val = met->evaluate_element(triPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //Check inverse mean ratio of ideal tri (INVERSE)
     val = imet->evaluate_element(triPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     
       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check mean ratio of ideal quad
     val = met->evaluate_element(quadPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //Check inverse mean ratio of ideal quad (INVERSE)
     val = imet->evaluate_element(quadPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check mean ratio of ideal tet
     val = met->evaluate_element(tetPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     //Check inverse mean ratio of ideal tet (INVERSE)
     val = imet->evaluate_element(tetPatch,&elems[0],err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check mean ratio of ideal hex
     val = met->evaluate_element(hexPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //Check inverse mean ratio of ideal hex (INVERSE)
     val = imet->evaluate_element(hexPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
   }

  void test_aspect_ratio_gamma()
   {
       //START WITH TRI's
     MsqError err;
     double val;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = AspectRatioGammaQualityMetric::create_new();
       //Check aspect ratio gamma of ideal tri
     val = met->evaluate_element(triPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     
       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check aspect ratio gamma of ideal tet
     val = met->evaluate_element(tetPatch,&elems[0],err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

   }
  
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "QualityMetricTest");

