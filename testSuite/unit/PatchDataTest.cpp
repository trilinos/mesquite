// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD: 19-Nov-02 at 15:12:42 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PatchDataTest.cpp

Unit testing of various functions in the PatchData class. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "PatchData.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;

class PatchDataTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(PatchDataTest);
   CPPUNIT_TEST (test_add_vertex);
   CPPUNIT_TEST_SUITE_END();
   
private:

   /* our 2D set up: 2 triangles and one quad are available
   ---------
   |\  |   |
   | \ |   |
   |  \|   |
   ---------
   */
   MsqVertex vtx_0_0;
   MsqVertex vtx_0_1;
   MsqVertex vtx_1_0;
   MsqVertex vtx_1_1;
   MsqVertex vtx_2_0;
   MsqVertex vtx_2_1;

   MsqMeshEntity tri1;
   MsqMeshEntity tri2;
   MsqMeshEntity quad1;   
   
   PatchData mPatch2D;

public:
  void setUp()
  {
     MsqError err;
     vtx_0_0.set(0,0,0);
     vtx_0_1.set(0,1,0);
     vtx_1_0.set(1,0,0);
     vtx_1_1.set(1,1,0);
     vtx_2_0.set(2,0,0);
     vtx_2_1.set(2,1,0);

     mPatch2D.reserve_vertex_capacity(6, err); MSQ_CHKERR(err);
     mPatch2D.add_vertex(NULL, NULL, 0,0,0, true, err); MSQ_CHKERR(err);
     mPatch2D.add_vertex(NULL, NULL, 0,1,0, true, err); MSQ_CHKERR(err);
     mPatch2D.add_vertex(NULL, NULL, 1,0,0, true, err); MSQ_CHKERR(err);
     mPatch2D.add_vertex(NULL, NULL, 1,1,0, true, err); MSQ_CHKERR(err);
     mPatch2D.add_vertex(NULL, NULL, 2,0,0, true, err); MSQ_CHKERR(err);
     mPatch2D.add_vertex(NULL, NULL, 2,1,0, true, err); MSQ_CHKERR(err);

     int ind[4];
     mPatch2D.reserve_element_capacity(3, err); MSQ_CHKERR(err);
     ind[0] = 0; ind[1]=2; ind[2]=1;
     mPatch2D.add_element(NULL, NULL, ind, TRIANGLE, err); MSQ_CHKERR(err);
     ind[0] = 1; ind[1]=2; ind[2]=3;
     mPatch2D.add_element(NULL, NULL, ind, TRIANGLE, err); MSQ_CHKERR(err);
     ind[0] = 3; ind[1]=2; ind[2]=4; ind[3]=5;
     mPatch2D.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);     
  }

  void tearDown()
  {
  }
  
public:
  PatchDataTest()
    {}
  
   void test_add_vertex()
   {
      MsqError err;
      
      MsqVertex* vertices = mPatch2D.get_vertex_array(err); MSQ_CHKERR(err);

      CPPUNIT_ASSERT_EQUAL(vtx_0_0, vertices[0]);
      CPPUNIT_ASSERT_EQUAL(vtx_0_1, vertices[1]);
      CPPUNIT_ASSERT_EQUAL(vtx_1_0, vertices[2]);
      CPPUNIT_ASSERT_EQUAL(vtx_1_1, vertices[3]);
      CPPUNIT_ASSERT_EQUAL(vtx_2_0, vertices[4]);
      CPPUNIT_ASSERT_EQUAL(vtx_2_1, vertices[5]);
   }

     
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "PatchDataTest");

