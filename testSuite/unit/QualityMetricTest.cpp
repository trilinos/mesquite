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
//  LAST-MOD: 13-Dec-02 at 08:26:36 by Thomas Leurent
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
#include "MultiplyQualityMetric.hpp"
#include "MsqMessage.hpp"
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
    //Test composite multiply
  CPPUNIT_TEST (test_composite_multiply);
    //Test averaging methods
  CPPUNIT_TEST (test_averaging_method);
  
  
  CPPUNIT_TEST_SUITE_END();
  
private:
  
  PatchData triPatch;
  PatchData quadPatch;
  PatchData tetPatch;
  PatchData hexPatch;
    //Tol used for double comparisons
  double qualTol;
  int pF;//PRINT_FLAG
public:
  void setUp()
  {
      //pF=1;//PRINT_FLAG IS ON
      pF=0;//PRINT_FLAG IS OFF
    qualTol = MSQ_MIN;
    MsqError err;
    size_t ind[20];
    
    size_t elem_ind[8];
    
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
     ind[8]=hexPatch.add_vertex(NULL, NULL, 2.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[9]=hexPatch.add_vertex(NULL, NULL, 2.0,1.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[10]=hexPatch.add_vertex(NULL, NULL, 2.0, -1.0, 1.0, true, err);
     MSQ_CHKERR(err);
     ind[11]=hexPatch.add_vertex(NULL, NULL, 3.0, 2.0, 1.0, true, err);
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
     bool v_flag;
       //START WITH TRI's
     MsqError err;
     double val, val2;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = ConditionNumberQualityMetric::create_new();
     ShapeQualityMetric *gmet = GeneralizedConditionNumberQualityMetric::create_new();
       //Check condition number of ideal tri
     met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tri.
     met->evaluate_element(triPatch,&elems[1],val,err); MSQ_CHKERR(err);
     gmet->evaluate_element(triPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN TRI %f", val2);
     
     CPPUNIT_ASSERT(fabs(val)<qualTol);
     
       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check condition number of ideal quad
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary quad.
     met->evaluate_element(quadPatch,&elems[1],val,err); MSQ_CHKERR(err);
     gmet->evaluate_element(quadPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN QUA %f", val2);
     
     CPPUNIT_ASSERT(fabs(val)<qualTol);

       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check condition number of ideal tet
     val = met->evaluate_element(tetPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tet.
     met->evaluate_element(tetPatch,&elems[1],val,err); MSQ_CHKERR(err);
     v_flag=gmet->evaluate_element(tetPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN TET %f", val2);
     
       //CPPUNIT_ASSERT(fabs(val)<qualTol);

       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check condition number of ideal hex
     met->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tet.
     v_flag=met->evaluate_element(hexPatch,&elems[1],val,err); MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nCON HEX %f", val);
     CPPUNIT_ASSERT(v_flag==true);
     
     gmet->evaluate_element(hexPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN HEX %f", val2);
     
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
     met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMEAN TRI %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //Check inverse mean ratio of ideal tri (INVERSE)
     imet->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nInv MEAN TRI %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     
       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check mean ratio of ideal quad
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMEAN QUAD %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //Check inverse mean ratio of ideal quad (INVERSE)
     imet->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nInv MEAN QUAD %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check mean ratio of ideal tet
     met->evaluate_element(tetPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMEAN TET %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     //Check inverse mean ratio of ideal tet (INVERSE)
     imet->evaluate_element(tetPatch,&elems[0],val,err); MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nInv MEAN TET %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check mean ratio of ideal hex
     met->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMEAN HEX %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
       //Check inverse mean ratio of ideal hex (INVERSE)
     imet->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nInv MEAN HEX %f", val);
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
     met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     
       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check aspect ratio gamma of ideal tet
     met->evaluate_element(tetPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

   }

  void test_composite_multiply()
   {
       //START WITH TRI's
     MsqError err;
     double val;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *mmet = MeanRatioQualityMetric::create_new();
     ShapeQualityMetric *cmet = ConditionNumberQualityMetric::create_new();
     CompositeQualityMetric *met = MultiplyQualityMetric::create_new(mmet,
                                                                     cmet,
                                                                     err);
       //Check ideal tri
     met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMULT TRI %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check ideal quad
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMULT QUAD %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check ideal tet
     met->evaluate_element(tetPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMULT TET %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);

       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check ideal hex
     met->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMULT HEX %f", val);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
   }
void test_averaging_method()
   {
       //USE QUAD's
     MsqError err;
     double val;
     MsqMeshEntity* elems;
     MsqVertex* verts = quadPatch.get_vertex_array(err);
     elems=quadPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = MeanRatioQualityMetric::create_new();
       //Check mean ratio of ideal quad
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     met->set_averaging_method(QualityMetric::GEOMETRIC, err);
       //Check mean ratio of ideal quad GEOMETRIC
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     met->set_averaging_method(QualityMetric::HARMONIC, err);
       //Check mean ratio of ideal quad HARMONIC
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     met->set_averaging_method(QualityMetric::LINEAR, err);
       //Check mean ratio of ideal quad LINEAR
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     met->set_averaging_method(QualityMetric::MAXIMUM, err);
       //Check mean ratio of ideal quad MAXIMUM
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     met->set_averaging_method(QualityMetric::MINIMUM, err);
       //Check mean ratio of ideal quad MINIMUM
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     met->set_averaging_method(QualityMetric::RMS, err);
       //Check mean ratio of ideal quad RMS
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-1.0)<qualTol);
     met->set_averaging_method(QualityMetric::SUM, err);
       //Check mean ratio of ideal SUM (NOTICE:: should be 4.0)
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(fabs(val-4.0)<qualTol);
     

   }

   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "QualityMetricTest");

