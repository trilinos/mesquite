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
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD:  9-Jun-04 at 14:50:51 by Thomas Leurent
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
#include "PatchDataInstances.hpp"

#include "cppunit/extensions/HelperMacros.h"

using namespace Mesquite;

using std::cout;
using std::cerr;
using std::endl;

class PatchDataTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PatchDataTest);
  CPPUNIT_TEST (test_num_corners);
  CPPUNIT_TEST (test_get_element_vertex_indices);
  CPPUNIT_TEST (test_get_vertex_element_indices);
  CPPUNIT_TEST (test_get_element_vertex_coordinates);
  CPPUNIT_TEST (test_move_free_vertices_constrained);
  CPPUNIT_TEST (test_movement_function);
  CPPUNIT_TEST (test_get_adj_elems_2d);
  CPPUNIT_TEST (test_get_minmax_element_area);
  CPPUNIT_TEST (test_get_barrier_delta);
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
     MsqPrintError err(cout);

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
     
     double coords[] = { 0,0,0,
                         0,1,0,
                         1,0,0,
                         1,1,0,
                         2,0,0,
                         2,1,0 };
     
     EntityTopology types[] = { TRIANGLE, TRIANGLE, QUADRILATERAL };
     
     size_t connectivity[] = { 0, 2, 1,
                               1, 2, 3,
                               3, 2, 4, 5 };
     
     size_t counts[] = { 3, 3, 4 };
     
     mPatch2D.fill( 6, coords,
                    3, types,
                    counts, connectivity,
                    0, err );
  }
  
  void tearDown()
  {
  }
  
public:
  PatchDataTest()
    {}
  
   void test_num_corners()
   {
     MsqPrintError err(cout);
     size_t n = mPatch2D.num_corners();
     CPPUNIT_ASSERT(n==10);
   }

   void test_get_element_vertex_indices()
   {

      MsqPrintError err(cout);
      
      std::vector<size_t> vtx_ind;
      std::vector<size_t> res;

      // test we get the right vertices for element 1 (tri)
      mPatch2D.get_element_vertex_indices(1, vtx_ind, err); CPPUNIT_ASSERT(!err);
      res.push_back(1); res.push_back(2); res.push_back(3);
      CPPUNIT_ASSERT( vtx_ind==res );

      // test we get the right vertices for element 2 (quad)
      vtx_ind.clear(); res.clear();
      mPatch2D.get_element_vertex_indices(2, vtx_ind, err); CPPUNIT_ASSERT(!err);
      res.push_back(3); res.push_back(2); res.push_back(4); res.push_back(5);
      CPPUNIT_ASSERT( vtx_ind==res );
   }

   void test_get_vertex_element_indices()
   {
     /*  1___3___5
         |\1|   |
         |0\| 2 |
         0---2---4   */
     MsqPrintError err(cout);
     
     std::vector<size_t> elem_ind;
     std::vector<size_t> res;
     
     mPatch2D.generate_vertex_to_element_data();
     
     // test we get the elements contiguous to vertex 3
     mPatch2D.get_vertex_element_indices(3, elem_ind,err); CPPUNIT_ASSERT(!err);
     res.push_back(2); res.push_back(1);
     CPPUNIT_ASSERT(res==elem_ind);
     
     // test we get the elements contiguous to vertex 2
     elem_ind.clear(); res.clear();
     mPatch2D.get_vertex_element_indices(2, elem_ind,err); CPPUNIT_ASSERT(!err);
     res.push_back(2); res.push_back(1); res.push_back(0);
     CPPUNIT_ASSERT(res==elem_ind);
   }

   void test_get_element_vertex_coordinates()
   {
      MsqPrintError err(cout);

      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(1, coords,err); CPPUNIT_ASSERT(!err);
      
      CPPUNIT_ASSERT( coords[0]==vtx_0_1 );
      CPPUNIT_ASSERT( coords[1]==vtx_1_0 );
      CPPUNIT_ASSERT( coords[2]==vtx_1_1 );
   }

   /* This tests the move_vertices() function as well as the
      PatchDataCoordsMemento functionality
      */
   void test_move_free_vertices_constrained()
   {
      MsqPrintError err(cout);

      // gets a memento of the patch coordinates.
      PatchDataVerticesMemento* coords_mem = mPatch2D.create_vertices_memento(err);
      CPPUNIT_ASSERT(!err);
      
      // Move the two first vertices in direction dk by step size s;
      Vector3D dk[6];
      dk[0].set(-1,-2,0);
      dk[1].set(-1, 2,0);
      double s = 0.3;
      mPatch2D.move_free_vertices_constrained(dk, 6, s, err); CPPUNIT_ASSERT(!err);

      // gets the new coordinates and  checks the vertices were displaced as expected.
      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      Vector3D new_vtx_0_0 = vtx_0_0 + s*dk[0];
      Vector3D new_vtx_0_1 = vtx_0_1 + s*dk[1];
      CPPUNIT_ASSERT(coords[0] == new_vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == new_vtx_0_1);

      // restore the PatchData to previous coords.
      mPatch2D.set_to_vertices_memento(coords_mem, err); CPPUNIT_ASSERT(!err);

      // gets the new coordinates and  checks the vertices are back to original.
      coords.clear();
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      CPPUNIT_ASSERT(coords[0] == vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == vtx_0_1);

      delete coords_mem;
   }

   void test_movement_function()
   {
      MsqPrintError err(cout);
      // gets a memento of the patch coordinates.
      PatchDataVerticesMemento* coords_mem = mPatch2D.create_vertices_memento(err);
      CPPUNIT_ASSERT(!err);
      
      // Move the two first vertices in direction dk by step size s;
      Vector3D dk[6];
      dk[0].set(0,-2,0);
      dk[1].set(-1,0,0);
      double s = 1;
      mPatch2D.move_free_vertices_constrained(dk, 6, 1, err); CPPUNIT_ASSERT(!err);
      // gets the new coordinates and  checks the vertices were displaced as expected.
      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      Vector3D new_vtx_0_0 = vtx_0_0 + s*dk[0];
      Vector3D new_vtx_0_1 = vtx_0_1 + s*dk[1];
      CPPUNIT_ASSERT(coords[0] == new_vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == new_vtx_0_1);
      double m_dist=mPatch2D.get_max_vertex_movement_squared(coords_mem,err);
      CPPUNIT_ASSERT(m_dist==4.0);
      // restore the PatchData to previous coords.
      mPatch2D.set_to_vertices_memento(coords_mem, err); CPPUNIT_ASSERT(!err);
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
   void test_get_adj_elems_2d()
   {
     MsqPrintError err(cout);
     std::vector<size_t> elems_0;
       //find elements sharing an edge with oth elem (should be 1)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 0, elems_0, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_0.back() == 1);
     std::vector<size_t> elems_1;
       //find elements sharing an edge with 1st elem (should be 0 and 2)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 1, elems_1, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_1.size() == 2);
     std::vector<size_t> elems_2;
       //find elements sharing an vert with 0th elem (should be 1 and 2).
     mPatch2D.get_adjacent_entities_via_n_dim(0, 0, elems_2, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_2.size() == 2);
     std::vector<size_t> elems_3;
     //find elements sharing an face with 0th elem (should be empty).
     mPatch2D.get_adjacent_entities_via_n_dim(2, 0, elems_3, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_3.size() == 0);
   }

  
   void test_get_minmax_element_area()
   {
     MsqPrintError err(cout);
     double min, max;
     mPatch2D.get_minmax_element_unsigned_area(min, max, err); CPPUNIT_ASSERT(!err);

     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, min, 0.0001 );
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, max, 0.0001 );
   }        
     
   void test_get_barrier_delta()
   {
     MsqPrintError err(cout);
     
     PatchData pd1;
     create_six_quads_patch_with_domain(pd1, err); CPPUNIT_ASSERT(!err);
     pd1.clear_computed_info();
     double b = pd1.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b, 0.00001 );
     destroy_patch_with_domain(pd1);

     PatchData pd2;
     create_six_quads_patch_inverted_with_domain(pd2, err); CPPUNIT_ASSERT(!err);
     pd2.clear_computed_info();
     b = pd2.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.003, b, 0.00001 );
     destroy_patch_with_domain(pd2);

     PatchData pd3;
     create_twelve_hex_patch(pd3, err); CPPUNIT_ASSERT(!err);
     pd3.clear_computed_info();
     b = pd3.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b, 0.00001 );

     PatchData pd4;
     create_twelve_hex_patch_inverted(pd4, err); CPPUNIT_ASSERT(!err);
     pd4.clear_computed_info();
     b = pd4.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0025, b, 0.00001 );
   } 
     
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "PatchDataTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "Unit");
