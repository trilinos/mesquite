/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file SamplePointsTest.cpp
 *  \brief Unit tests for SamplePoints class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "SamplePoints.hpp"
#include "TopologyInfo.hpp"
#include <cppunit/extensions/HelperMacros.h>

using namespace Mesquite;

class SamplePointsTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(SamplePointsTest);
  CPPUNIT_TEST (test_get_set_all);
  CPPUNIT_TEST (test_initialization);
  CPPUNIT_TEST (test_get_set_clear);
  CPPUNIT_TEST_SUITE_END();
  
public:
  
  void test_get_set_all ();
  void test_initialization ();
  void test_get_set_clear ();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SamplePointsTest, "SamplePointsTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SamplePointsTest, "Unit");

#define TEST_GET_SET_TYPE( T ) \
  def = sp.get_sample_topology_bits( T ); \
  sp.set_sample_topology_bits( T, 0u ); \
  cur = sp.get_sample_topology_bits( T ); \
  CPPUNIT_ASSERT_EQUAL( 0u, cur ); \
  sp.set_sample_topology_bits( T, def ); \
  cur = sp.get_sample_topology_bits( T ); \
  CPPUNIT_ASSERT_EQUAL( def, cur )

void SamplePointsTest::test_get_set_all()
{
  SamplePoints sp( true );
  unsigned def, cur;
  
  TEST_GET_SET_TYPE( TRIANGLE );
  TEST_GET_SET_TYPE( QUADRILATERAL );
  TEST_GET_SET_TYPE( PRISM );
  TEST_GET_SET_TYPE( PYRAMID );
}

#define TEST_INIT_BITS( corners, edges, faces, elem, type ) \
  do { \
    SamplePoints def_sp( corners, edges, faces, elem ); \
    CPPUNIT_ASSERT_EQUAL( corners, def_sp.will_sample_at( type, 0 ) ); \
    CPPUNIT_ASSERT_EQUAL(   edges, def_sp.will_sample_at( type, 1 ) ); \
    CPPUNIT_ASSERT_EQUAL(   faces, def_sp.will_sample_at( type, 2 ) ); \
    CPPUNIT_ASSERT_EQUAL(    elem, def_sp.will_sample_at( type, 3 ) ); \
  } while (false)

void SamplePointsTest::test_initialization()
{
  SamplePoints def(true), none(false);
  
  /* Test that all element types have a default value */
  CPPUNIT_ASSERT( def.get_sample_topology_bits( TRIANGLE ) != 0 );
  CPPUNIT_ASSERT( def.get_sample_topology_bits( QUADRILATERAL ) != 0 );
  CPPUNIT_ASSERT( def.get_sample_topology_bits( TETRAHEDRON ) != 0 );
  CPPUNIT_ASSERT( def.get_sample_topology_bits( PYRAMID ) != 0 );
  CPPUNIT_ASSERT( def.get_sample_topology_bits( PRISM ) != 0 );
  CPPUNIT_ASSERT( def.get_sample_topology_bits( HEXAHEDRON ) != 0 );
  
  /* Test that no bits are set if defaults were not requested */
  CPPUNIT_ASSERT_EQUAL( 0u, none.get_sample_topology_bits( TRIANGLE ) );
  CPPUNIT_ASSERT_EQUAL( 0u, none.get_sample_topology_bits( QUADRILATERAL ) );
  CPPUNIT_ASSERT_EQUAL( 0u, none.get_sample_topology_bits( TETRAHEDRON ) );
  CPPUNIT_ASSERT_EQUAL( 0u, none.get_sample_topology_bits( PYRAMID ) );
  CPPUNIT_ASSERT_EQUAL( 0u, none.get_sample_topology_bits( PRISM ) );
  CPPUNIT_ASSERT_EQUAL( 0u, none.get_sample_topology_bits( HEXAHEDRON ) );
  
  /* Test clear_all() method */
  def.clear_all();
  CPPUNIT_ASSERT_EQUAL( 0u, def.get_sample_topology_bits( TRIANGLE ) );
  CPPUNIT_ASSERT_EQUAL( 0u, def.get_sample_topology_bits( QUADRILATERAL ) );
  CPPUNIT_ASSERT_EQUAL( 0u, def.get_sample_topology_bits( TETRAHEDRON ) );
  CPPUNIT_ASSERT_EQUAL( 0u, def.get_sample_topology_bits( PYRAMID ) );
  CPPUNIT_ASSERT_EQUAL( 0u, def.get_sample_topology_bits( PRISM ) );
  CPPUNIT_ASSERT_EQUAL( 0u, def.get_sample_topology_bits( HEXAHEDRON ) );
  
  TEST_INIT_BITS( true, false, false, false, TRIANGLE );
  TEST_INIT_BITS( true, false, false, false, HEXAHEDRON );
  TEST_INIT_BITS( false, true, false, false, TRIANGLE );
  TEST_INIT_BITS( false, true, false, false, HEXAHEDRON );
  TEST_INIT_BITS( false, false, true, false, TRIANGLE );
  TEST_INIT_BITS( false, false, true, false, HEXAHEDRON );
  TEST_INIT_BITS( false, false, false, true, HEXAHEDRON );
  TEST_INIT_BITS( true, true, false, false, TRIANGLE );
  TEST_INIT_BITS( true, true, false, false, HEXAHEDRON );
  TEST_INIT_BITS( true, false, true, false, TRIANGLE );
  TEST_INIT_BITS( true, false, true, false, HEXAHEDRON );
  TEST_INIT_BITS( true, false, false, true, HEXAHEDRON );
}

static void test_get_set_clear_dim( unsigned dim, EntityTopology type )
{
  SamplePoints sp(false);
  CPPUNIT_ASSERT_EQUAL( false, sp.will_sample_at( type, dim ) );
  sp.sample_at( type, dim );
  CPPUNIT_ASSERT_EQUAL( true, sp.will_sample_at( type, dim ) );
  CPPUNIT_ASSERT_EQUAL( false, sp.will_sample_at( type, (dim ? 0 : 1) ) );
  sp.dont_sample_at( type, dim );
  CPPUNIT_ASSERT_EQUAL( false, sp.will_sample_at( type, dim ) );
}  

void SamplePointsTest::test_get_set_clear()
{
  SamplePoints sp(false);
  test_get_set_clear_dim( 0, TRIANGLE );
  test_get_set_clear_dim( 1, TRIANGLE );
  test_get_set_clear_dim( 2, TRIANGLE );
  test_get_set_clear_dim( 0, PYRAMID );
  test_get_set_clear_dim( 1, PYRAMID );
  test_get_set_clear_dim( 2, PYRAMID );
  test_get_set_clear_dim( 3, PYRAMID );
}

