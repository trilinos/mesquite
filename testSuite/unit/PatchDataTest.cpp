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
//  LAST-MOD: 20-Jan-03 at 16:21:12 by Thomas Leurent
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
   CPPUNIT_TEST (test_get_element_vertex_coordinates);
   CPPUNIT_TEST (test_move_vertices);
   CPPUNIT_TEST (test_get_adj_elems_2d);
  
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

     /* our 2D set up: 2 triangles and one quad are available
       1___3___5
        |\1|   |
        |0\| 2 |
       0---2---4
     */
     vtx_0_0.set(0,0,0);
     vtx_0_1.set(0,1,0);
     vtx_1_0.set(1,0,0);
     vtx_1_1.set(1,1,0);
     vtx_2_0.set(2,0,0);
     vtx_2_1.set(2,1,0);

     mPatch2D.reserve_vertex_capacity(6, err); MSQ_CHKERR(err);
     mPatch2D.add_vertex(NULL, NULL, 0,0,0, true, err); 
     mPatch2D.add_vertex(NULL, NULL, 0,1,0, true, err); 
     mPatch2D.add_vertex(NULL, NULL, 1,0,0, true, err, MsqVertex::MSQ_HARD_FIXED); 
     mPatch2D.add_vertex(NULL, NULL, 1,1,0, true, err, MsqVertex::MSQ_HARD_FIXED); 
     mPatch2D.add_vertex(NULL, NULL, 2,0,0, true, err, MsqVertex::MSQ_HARD_FIXED); 
     mPatch2D.add_vertex(NULL, NULL, 2,1,0, true, err, MsqVertex::MSQ_HARD_FIXED); 
     MSQ_CHKERR(err);

     size_t ind[4];
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
  
#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_reserve_vertex_capacity" 
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
   
#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_reserve_element_capacity" 
   void test_reserve_element_capacity()
   {
      MsqError err;

      size_t ind[3];
      ind[0] = 4; ind[1]=2; ind[2]=1;
      mPatch2D.add_element(NULL, NULL, ind, TRIANGLE, err);
      CPPUNIT_ASSERT(err.errorOn);
      err.reset();
      
      int num = mPatch2D.num_elements();
      CPPUNIT_ASSERT_EQUAL(3, num);
   }
   
#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_add_vertex" 
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

      size_t ind = mPatch2D.add_vertex(NULL, NULL, 0,0,0, true, err);
      size_t zero = 0;
      CPPUNIT_ASSERT_EQUAL(zero, ind);
      CPPUNIT_ASSERT_MESSAGE("vertex already existed.", !err.errorOn);
   }

#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_get_element_vertex_indices" 
   void test_get_element_vertex_indices()
   {

      MsqError err;
      
      std::vector<size_t> vtx_ind;
      std::vector<size_t> res;

      // test we get the right vertices for element 1 (tri)
      mPatch2D.get_element_vertex_indices(1, vtx_ind, err); MSQ_CHKERR(err);
      res.push_back(1); res.push_back(2); res.push_back(3);
      CPPUNIT_ASSERT( vtx_ind==res );

      // test we get the right vertices for element 2 (quad)
      vtx_ind.clear(); res.clear();
      mPatch2D.get_element_vertex_indices(2, vtx_ind, err); MSQ_CHKERR(err);
      res.push_back(3); res.push_back(2); res.push_back(4); res.push_back(5);
      CPPUNIT_ASSERT( vtx_ind==res );
   }

#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_get_vertex_element_indices" 
   void test_get_vertex_element_indices()
   {
     /*  1___3___5
         |\1|   |
         |0\| 2 |
         0---2---4   */
     MsqError err;
     
     std::vector<size_t> elem_ind;
     std::vector<size_t> res;
     
     mPatch2D.generate_vertex_to_element_data();
     
     // test we get the elements contiguous to vertex 3
     mPatch2D.get_vertex_element_indices(3, elem_ind,err); MSQ_CHKERR(err);
     res.push_back(2); res.push_back(1);
     CPPUNIT_ASSERT(res==elem_ind);
     
     // test we get the elements contiguous to vertex 2
     elem_ind.clear(); res.clear();
     mPatch2D.get_vertex_element_indices(2, elem_ind,err); MSQ_CHKERR(err);
     res.push_back(2); res.push_back(1); res.push_back(0);
     CPPUNIT_ASSERT(res==elem_ind);
   }

#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_get_element_vertex_coordinates" 
   void test_get_element_vertex_coordinates()
   {
      MsqError err;

      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(1, coords,err); MSQ_CHKERR(err);
      
      CPPUNIT_ASSERT( coords[0]==vtx_0_1 );
      CPPUNIT_ASSERT( coords[1]==vtx_1_0 );
      CPPUNIT_ASSERT( coords[2]==vtx_1_1 );
   }

   /* This tests the move_vertices() function as well as the
      PatchDataCoordsMemento functionality
      */
#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_move_vertices" 
   void test_move_vertices()
   {
      MsqError err;

      // gets a memento of the patch coordinates.
      PatchDataVerticesMemento* coords_mem = mPatch2D.create_vertices_memento(err);
      MSQ_CHKERR(err);
      
      // Move the two first vertices in direction dk by step size s;
      Vector3D dk[2];
      dk[0].set(-1,-2,0);
      dk[1].set(-1, 2,0);
      double s = 0.3;
      mPatch2D.move_vertices(dk, 2, s, err); MSQ_CHKERR(err);

      // gets the new coordinates and  checks the vertices were displaced as expected.
      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      Vector3D new_vtx_0_0 = vtx_0_0 + s*dk[0];
      Vector3D new_vtx_0_1 = vtx_0_1 + s*dk[1];
      CPPUNIT_ASSERT(coords[0] == new_vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == new_vtx_0_1);

      // restore the PatchData to previous coords.
      mPatch2D.set_to_vertices_memento(coords_mem, err); MSQ_CHKERR(err);

      // gets the new coordinates and  checks the vertices are back to original.
      coords.clear();
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      CPPUNIT_ASSERT(coords[0] == vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == vtx_0_1);

      delete coords_mem;
   }
/*Tests the function PatchData::get_adjacent_entities_via_n_dim()
  which finds the elements adjacent to a given element.  If 'n'
  equals 0 the elements must share a vertex; if 'n' equals 1 the
  elements must share an edge; and if 'n' equals 2 the elements
  must share a face.*/
#undef __FUNC__
#define __FUNC__ "PatchDataTest::test_get_adj_elems_2d" 
   void test_get_adj_elems_2d()
   {
     MsqError err;
     std::vector<size_t> elems_0;
       //find elements sharing an edge with oth elem (should be 1)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 0, elems_0, err);
     MSQ_CHKERR(err);
     CPPUNIT_ASSERT(elems_0.back() == 1);
     std::vector<size_t> elems_1;
       //find elements sharing an edge with 1st elem (should be 0 and 2)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 1, elems_1, err);
     MSQ_CHKERR(err);
     CPPUNIT_ASSERT(elems_1.size() == 2);
     std::vector<size_t> elems_2;
       //find elements sharing an vert with 0th elem (should be 1 and 2).
     mPatch2D.get_adjacent_entities_via_n_dim(0, 0, elems_2, err);
     MSQ_CHKERR(err);
     CPPUNIT_ASSERT(elems_2.size() == 2);
     std::vector<size_t> elems_3;
     //find elements sharing an face with 0th elem (should be empty).
     mPatch2D.get_adjacent_entities_via_n_dim(2, 0, elems_3, err);
     MSQ_CHKERR(err);
     CPPUNIT_ASSERT(elems_3.size() == 0);
   }
  
     
        
     
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "PatchDataTest");

