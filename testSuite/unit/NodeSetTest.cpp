/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file NodeSetTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "NodeSet.hpp"
#include "UnitUtil.hpp"

using namespace Mesquite;
using namespace std;

class NodeSetTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE( NodeSetTest );
    CPPUNIT_TEST(test_init);
    CPPUNIT_TEST(test_clear);
    CPPUNIT_TEST(test_set_get_clear_simple);
    CPPUNIT_TEST(test_set_get_clear_dim);
    CPPUNIT_TEST(test_set_node);
    CPPUNIT_TEST(test_clear_node);
    CPPUNIT_TEST(test_num_nodes);
    CPPUNIT_TEST(test_have_any);
    CPPUNIT_TEST(test_set_all);
    CPPUNIT_TEST(test_clear_all);
    CPPUNIT_TEST(test_num_before);
    CPPUNIT_TEST_SUITE_END();
    
  public:

    void test_init();
    void test_clear();
    void test_set_get_clear_simple();
    void test_set_get_clear_dim();
    void test_set_node();
    void test_clear_node();
    void test_num_nodes();
    void test_have_any();
    void test_set_all();
    void test_clear_all();
    void test_num_before();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(NodeSetTest, "NodeSetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(NodeSetTest, "Unit");

void NodeSetTest::test_init()
{
  NodeSet set;
  CPPUNIT_ASSERT( !set.get_bits() );
}

void NodeSetTest::test_clear()
{
  NodeSet set;
  set.set_all_corner_nodes();
  set.clear();
  CPPUNIT_ASSERT( !set.get_bits() );
  set.set_all_mid_edge_nodes();
  set.clear();
  CPPUNIT_ASSERT( !set.get_bits() );
  set.set_all_mid_nodes();
  set.clear();
  CPPUNIT_ASSERT( !set.get_bits() );
}

void NodeSetTest::test_set_get_clear_simple()
{
  NodeSet set;
  for (unsigned i = 0; i < NodeSet::NUM_CORNER_BITS; ++i) {
    CPPUNIT_ASSERT( !set.corner_node(i) );
    set.set_corner_node( i );
    CPPUNIT_ASSERT( set.corner_node(i) );
    set.clear_corner_node( i );
    CPPUNIT_ASSERT( !set.corner_node(i) );
  }
  for (unsigned i = 0; i < NodeSet::NUM_EDGE_BITS; ++i) {
    CPPUNIT_ASSERT( !set.mid_edge_node(i) );
    set.set_mid_edge_node( i );
    CPPUNIT_ASSERT( set.mid_edge_node(i) );
    set.clear_mid_edge_node( i );
    CPPUNIT_ASSERT( !set.mid_edge_node(i) );
  }
  for (unsigned i = 0; i < NodeSet::NUM_FACE_BITS; ++i) {
    CPPUNIT_ASSERT( !set.mid_face_node(i) );
    set.set_mid_face_node( i );
    CPPUNIT_ASSERT( set.mid_face_node(i) );
    set.clear_mid_face_node( i );
    CPPUNIT_ASSERT( !set.mid_face_node(i) );
  }
  for (unsigned i = 0; i < NodeSet::NUM_REGION_BITS; ++i) {
    CPPUNIT_ASSERT( !set.mid_region_node(i) );
    set.set_mid_region_node( i );
    CPPUNIT_ASSERT( set.mid_region_node(i) );
    set.clear_mid_region_node( i );
    CPPUNIT_ASSERT( !set.mid_region_node(i) );
  }
}

void NodeSetTest::test_set_get_clear_dim()
{
  NodeSet set;
  for (unsigned i = 0; i < NodeSet::NUM_CORNER_BITS; ++i) {
    CPPUNIT_ASSERT( !set.node(Sample(0,i)) );
    set.set_node( Sample(0, i) );
    CPPUNIT_ASSERT( set.node(Sample(0,i)) );
    set.clear_node( Sample(0, i) );
    CPPUNIT_ASSERT( !set.node(Sample(0,i)) );
  }
  for (unsigned i = 0; i < NodeSet::NUM_EDGE_BITS; ++i) {
    CPPUNIT_ASSERT( !set.node(Sample(1,i)) );
    set.set_node( Sample(1, i) );
    CPPUNIT_ASSERT( set.node(Sample(1,i)) );
    set.clear_node( Sample(1, i) );
    CPPUNIT_ASSERT( !set.node(Sample(1,i)) );
  }
  for (unsigned i = 0; i < NodeSet::NUM_FACE_BITS; ++i) {
    CPPUNIT_ASSERT( !set.node(Sample(2,i)) );
    set.set_node( Sample(2, i) );
    CPPUNIT_ASSERT( set.node(Sample(2,i)) );
    set.clear_node( Sample(2, i) );
    CPPUNIT_ASSERT( !set.node(Sample(2,i)) );
  }
  for (unsigned i = 0; i < NodeSet::NUM_REGION_BITS; ++i) {
    CPPUNIT_ASSERT( !set.node(Sample(3,i)) );
    set.set_node( Sample(3, i) );
    CPPUNIT_ASSERT( set.node(Sample(3,i)) );
    set.clear_node( Sample(3, i) );
    CPPUNIT_ASSERT( !set.node(Sample(3,i)) );
  }
}

void NodeSetTest::test_set_node()
{
  NodeSet set;
  for (unsigned i = 0; i < NodeSet::NUM_CORNER_BITS; ++i) {
    set.set_corner_node( i );
    CPPUNIT_ASSERT_EQUAL( 1u << (NodeSet::CORNER_OFFSET + i), set.get_bits() );
    set.clear();
  }
  for (unsigned i = 0; i < NodeSet::NUM_EDGE_BITS; ++i) {
    set.set_mid_edge_node( i );
    CPPUNIT_ASSERT_EQUAL( 1u << (NodeSet::EDGE_OFFSET + i), set.get_bits() );
    set.clear();
  }
  for (unsigned i = 0; i < NodeSet::NUM_FACE_BITS; ++i) {
    set.set_mid_face_node( i );
    CPPUNIT_ASSERT_EQUAL( 1u << (NodeSet::FACE_OFFSET + i), set.get_bits() );
    set.clear();
  }
  for (unsigned i = 0; i < NodeSet::NUM_REGION_BITS; ++i) {
    set.set_mid_region_node( i );
    CPPUNIT_ASSERT_EQUAL( 1u << (NodeSet::REGION_OFFSET + i), set.get_bits() );
    set.clear();
  }
}

void NodeSetTest::test_clear_node()
{
  NodeSet set;
  for (unsigned i = 0; i < NodeSet::NUM_CORNER_BITS; ++i) {
    set.set_all_nodes();
    set.clear_corner_node( i );
    CPPUNIT_ASSERT_EQUAL( ~(1u << (NodeSet::CORNER_OFFSET + i)), set.get_bits() );
  }
  for (unsigned i = 0; i < NodeSet::NUM_EDGE_BITS; ++i) {
    set.set_all_nodes();
    set.clear_mid_edge_node( i );
    CPPUNIT_ASSERT_EQUAL( ~(1u << (NodeSet::EDGE_OFFSET + i)), set.get_bits() );
  }
  for (unsigned i = 0; i < NodeSet::NUM_FACE_BITS; ++i) {
    set.set_all_nodes();
    set.clear_mid_face_node( i );
    CPPUNIT_ASSERT_EQUAL( ~(1u << (NodeSet::FACE_OFFSET + i)), set.get_bits() );
  }
  for (unsigned i = 0; i < NodeSet::NUM_REGION_BITS; ++i) {
    set.set_all_nodes();
    set.clear_mid_region_node( i );
    CPPUNIT_ASSERT_EQUAL( ~(1u << (NodeSet::REGION_OFFSET + i)), set.get_bits() );
  }
}

void NodeSetTest::test_num_nodes()
{
  NodeSet set;
  set.clear();
  set.set_corner_node(1);
  CPPUNIT_ASSERT_EQUAL( 1u, set.num_nodes() );
  set.set_corner_node(3);
  CPPUNIT_ASSERT_EQUAL( 2u, set.num_nodes() );
  set.set_mid_region_node(0);
  CPPUNIT_ASSERT_EQUAL( 3u, set.num_nodes() );
  set.set_mid_edge_node(1);
  set.set_mid_edge_node(2);
  CPPUNIT_ASSERT_EQUAL( 5u, set.num_nodes() );
  set.set_all_nodes();
  CPPUNIT_ASSERT_EQUAL( (unsigned)NodeSet::NUM_TOTAL_BITS, set.num_nodes() );
}

void NodeSetTest::test_have_any()
{
  NodeSet set;
  set.clear();
  
  CPPUNIT_ASSERT( !set.have_any_corner_node() );
  set.set_corner_node(0);
  CPPUNIT_ASSERT( set.have_any_corner_node() );
  set.clear_corner_node(0);
  CPPUNIT_ASSERT( !set.have_any_corner_node() );
  set.set_corner_node( NodeSet::NUM_CORNER_BITS-1 );
  CPPUNIT_ASSERT( set.have_any_corner_node() );
  CPPUNIT_ASSERT( !set.have_any_mid_node() );

  CPPUNIT_ASSERT( !set.have_any_mid_edge_node() );
  set.set_mid_edge_node(0);
  CPPUNIT_ASSERT( set.have_any_mid_edge_node() );
  CPPUNIT_ASSERT( set.have_any_mid_node() );
  set.clear_mid_edge_node(0);
  CPPUNIT_ASSERT( !set.have_any_mid_edge_node() );
  CPPUNIT_ASSERT( !set.have_any_mid_node() );
  set.set_mid_edge_node( NodeSet::NUM_EDGE_BITS-1 );
  CPPUNIT_ASSERT( set.have_any_mid_edge_node() );
  CPPUNIT_ASSERT( set.have_any_mid_node() );
  set.clear();

  CPPUNIT_ASSERT( !set.have_any_mid_face_node() );
  set.set_mid_face_node(0);
  CPPUNIT_ASSERT( set.have_any_mid_face_node() );
  CPPUNIT_ASSERT( set.have_any_mid_node() );
  set.clear_mid_face_node(0);
  CPPUNIT_ASSERT( !set.have_any_mid_face_node() );
  CPPUNIT_ASSERT( !set.have_any_mid_node() );
  set.set_mid_face_node( NodeSet::NUM_FACE_BITS-1 );
  CPPUNIT_ASSERT( set.have_any_mid_face_node() );
  CPPUNIT_ASSERT( set.have_any_mid_node() );
  set.clear();

  CPPUNIT_ASSERT( !set.have_any_mid_region_node() );
  set.set_mid_region_node(0);
  CPPUNIT_ASSERT( set.have_any_mid_region_node() );
  CPPUNIT_ASSERT( set.have_any_mid_node() );
  set.clear_mid_region_node(0);
  CPPUNIT_ASSERT( !set.have_any_mid_region_node() );
  CPPUNIT_ASSERT( !set.have_any_mid_node() );
  set.set_mid_region_node( NodeSet::NUM_REGION_BITS-1 );
  CPPUNIT_ASSERT( set.have_any_mid_region_node() );
  CPPUNIT_ASSERT( set.have_any_mid_node() );
  set.clear();
}

static bool check_all( NodeSet set, unsigned dim, bool value )
{
  unsigned count = 0;
  switch (dim) { case 0: count = NodeSet::NUM_CORNER_BITS; break;
                 case 1: count = NodeSet::NUM_EDGE_BITS; break;
                 case 2: count = NodeSet::NUM_FACE_BITS; break;
                 case 3: count = NodeSet::NUM_REGION_BITS; break;
                 default: CPPUNIT_ASSERT(false); }
  for (unsigned i = 0; i < count; ++i)
    if (set.node( Sample(dim, i) ) != value )
      return false;
  return true;
}
  
void NodeSetTest::test_set_all()
{
  NodeSet set;
  set.set_all_corner_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, true ) );
  CPPUNIT_ASSERT( check_all( set, 1, false ) );
  CPPUNIT_ASSERT( check_all( set, 2, false ) );
  CPPUNIT_ASSERT( check_all( set, 3, false ) );
  set.clear();
  set.set_all_mid_edge_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, false ) );
  CPPUNIT_ASSERT( check_all( set, 1, true ) );
  CPPUNIT_ASSERT( check_all( set, 2, false ) );
  CPPUNIT_ASSERT( check_all( set, 3, false ) );
  set.clear();
  set.set_all_mid_face_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, false ) );
  CPPUNIT_ASSERT( check_all( set, 1, false ) );
  CPPUNIT_ASSERT( check_all( set, 2, true ) );
  CPPUNIT_ASSERT( check_all( set, 3, false ) );
  set.clear();
  set.set_all_mid_region_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, false ) );
  CPPUNIT_ASSERT( check_all( set, 1, false ) );
  CPPUNIT_ASSERT( check_all( set, 2, false ) );
  CPPUNIT_ASSERT( check_all( set, 3, true ) );
}

void NodeSetTest::test_clear_all()
{
  NodeSet set;
  set.set_all_nodes();
  set.clear_all_corner_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, false ) );
  CPPUNIT_ASSERT( check_all( set, 1, true ) );
  CPPUNIT_ASSERT( check_all( set, 2, true ) );
  CPPUNIT_ASSERT( check_all( set, 3, true ) );
  set.set_all_nodes();
  set.clear_all_mid_edge_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, true ) );
  CPPUNIT_ASSERT( check_all( set, 1, false ) );
  CPPUNIT_ASSERT( check_all( set, 2, true ) );
  CPPUNIT_ASSERT( check_all( set, 3, true ) );
  set.set_all_nodes();
  set.clear_all_mid_face_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, true ) );
  CPPUNIT_ASSERT( check_all( set, 1, true ) );
  CPPUNIT_ASSERT( check_all( set, 2, false ) );
  CPPUNIT_ASSERT( check_all( set, 3, true ) );
  set.set_all_nodes();
  set.clear_all_mid_region_nodes();
  CPPUNIT_ASSERT( check_all( set, 0, true ) );
  CPPUNIT_ASSERT( check_all( set, 1, true ) );
  CPPUNIT_ASSERT( check_all( set, 2, true ) );
  CPPUNIT_ASSERT( check_all( set, 3, false ) );
}

void NodeSetTest::test_num_before()
{
  NodeSet set;
  set.clear();
  set.set_mid_face_node( 2 );
  CPPUNIT_ASSERT_EQUAL( 0u, set.num_before_mid_face( 2 ) );
  CPPUNIT_ASSERT_EQUAL( 1u, set.num_before_mid_face( 3 ) );
  set.set_corner_node( 0 );
  CPPUNIT_ASSERT_EQUAL( 1u, set.num_before_mid_face( 2 ) );
  CPPUNIT_ASSERT_EQUAL( 2u, set.num_before_mid_face( 3 ) );
  CPPUNIT_ASSERT_EQUAL( 0u, set.num_before_corner( 0 ) );
  CPPUNIT_ASSERT_EQUAL( 1u, set.num_before_corner( 1 ) );
  set.set_all_corner_nodes();
  CPPUNIT_ASSERT_EQUAL( (unsigned)NodeSet::NUM_CORNER_BITS, set.num_before_mid_edge(0) );
  CPPUNIT_ASSERT_EQUAL( (unsigned)NodeSet::NUM_CORNER_BITS-1, set.num_before_corner(NodeSet::NUM_CORNER_BITS-1) );
  CPPUNIT_ASSERT_EQUAL( 0u, set.num_before_corner(0) );
  CPPUNIT_ASSERT_EQUAL( 1u, set.num_before_corner(1) );
  CPPUNIT_ASSERT_EQUAL( 2u, set.num_before_corner(2) );
  CPPUNIT_ASSERT_EQUAL( 3u, set.num_before_corner(3) );
  CPPUNIT_ASSERT_EQUAL( 4u, set.num_before_corner(4) );
  set.set_all_nodes();
  CPPUNIT_ASSERT_EQUAL( (unsigned)NodeSet::NUM_TOTAL_BITS - 1, set.num_before_mid_region( NodeSet::NUM_REGION_BITS -1 ) );
  set.clear_mid_edge_node( 0 );
  CPPUNIT_ASSERT_EQUAL( (unsigned)NodeSet::NUM_TOTAL_BITS - 2, set.num_before_mid_region( NodeSet::NUM_REGION_BITS -1 ) );
  set.clear_mid_edge_node( 1 );
  CPPUNIT_ASSERT_EQUAL( (unsigned)NodeSet::NUM_TOTAL_BITS - 3, set.num_before_mid_region( NodeSet::NUM_REGION_BITS -1 ) );
}
