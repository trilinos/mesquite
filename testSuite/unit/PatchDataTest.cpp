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
//  LAST-MOD: 20-Nov-02 at 10:52:57 by Thomas Leurent
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
   CPPUNIT_TEST (test_reserve_vertex_capacity);
   CPPUNIT_TEST (test_reserve_element_capacity);
   CPPUNIT_TEST (test_get_element_vertex_indices);
   CPPUNIT_TEST (test_get_vertex_element_indices);   
   CPPUNIT_TEST_SUITE_END();
   
private:

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

     /* our 2D set up: 2 triangles and one quad are available
        
       1___3___5
        |\2|   |
        |1\| 3 |
       0---2---4
     */
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
  
   void test_reserve_vertex_capacity()
   {
      MsqError err;

      // checks that adding too many vertices returns an error.
      mPatch2D.add_vertex(NULL, NULL, 1.5,5.6,0, true, err);
      CPPUNIT_ASSERT(err.errorOn);
      err.reset();

      mPatch2D.add_vertex(NULL, NULL, 1.5,5.6,0, false, err);
      CPPUNIT_ASSERT(err.errorOn);
      err.reset();

      int num = mPatch2D.num_vertices();
      CPPUNIT_ASSERT_EQUAL(6, num);
   }
   
   void test_reserve_element_capacity()
   {
      MsqError err;

      int ind[3];
      ind[0] = 4; ind[1]=2; ind[2]=1;
      mPatch2D.add_element(NULL, NULL, ind, TRIANGLE, err);
      CPPUNIT_ASSERT(err.errorOn);
      err.reset();
      
      int num = mPatch2D.num_elements();
      CPPUNIT_ASSERT_EQUAL(3, num);
   }
   
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

      int ind = mPatch2D.add_vertex(NULL, NULL, 0,0,0, true, err);
      CPPUNIT_ASSERT_EQUAL(0, ind);
      CPPUNIT_ASSERT_MESSAGE("vertex already existed.", !err.errorOn);
   }

   void test_get_element_vertex_indices() {

      MsqError err;
      
      std::vector<size_t> vtx_ind;
      std::vector<size_t> res;
      
      mPatch2D.get_element_vertex_indices(1, vtx_ind, err); MSQ_CHKERR(err);
      res.push_back(1); res.push_back(2); res.push_back(3);
      CPPUNIT_ASSERT( vtx_ind==res );

      vtx_ind.clear(); res.clear();
      mPatch2D.get_element_vertex_indices(2, vtx_ind, err); MSQ_CHKERR(err);
      res.push_back(3); res.push_back(2); res.push_back(4); res.push_back(5);
      CPPUNIT_ASSERT( vtx_ind==res );
   }

   void test_get_vertex_element_indices() {
      MsqError err;
      
      std::vector<size_t> elem_ind;
      std::vector<size_t> res;
      
      mPatch2D.get_vertex_element_indices(3, elem_ind);
      res.push_back(2); res.push_back(3);
      CPPUNIT_ASSERT(res==elem_ind);
      
      elem_ind.clear(); res.clear();
      mPatch2D.get_vertex_element_indices(2, elem_ind);
      res.push_back(1); res.push_back(2); res.push_back(3);
      CPPUNIT_ASSERT(res==elem_ind);
   }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "PatchDataTest");

