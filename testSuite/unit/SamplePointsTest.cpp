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
  CPPUNIT_TEST (test_corners);
  CPPUNIT_TEST (test_edges);
  CPPUNIT_TEST (test_corners_and_edges);
  CPPUNIT_TEST (test_edges_and_faces);
  CPPUNIT_TEST (test_elem );
  CPPUNIT_TEST (test_corners_and_elem );
  CPPUNIT_TEST_SUITE_END();
  
public:
  
  void test_get_set_all ();
  void test_initialization ();
  void test_get_set_clear ();
  void test_corners ();
  void test_edges ();
  void test_corners_and_edges ();
  void test_edges_and_faces ();
  void test_elem ();
  void test_corners_and_elem ();
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

static void test_topo_common( bool corners, 
                              bool edges,
                              bool faces,
                              bool elem,
                              EntityTopology type )
{
    // expected number of samples of each dimension
  unsigned expected_by_dim[4] = {0, 0, 0, 0};
  unsigned seen_by_dim[4] = { 0, 0, 0, 0 };
    // bool for each possible sample, indexed by dmension
  std::vector<bool> seen[4];
  
    // initialize everything
  seen[0].resize( TopologyInfo::corners(type), false );
  seen[1].resize( TopologyInfo::edges(type), false );
  if (corners) 
    expected_by_dim[0] = type == PYRAMID ? 4 : TopologyInfo::corners( type );
  if (edges)
    expected_by_dim[1] = TopologyInfo::edges( type );
  if (TopologyInfo::dimension(type) == 2) {
    expected_by_dim[2] = faces ? 1 : 0;
    seen[2].resize( 1, false );
  }
  else {
    if (faces)
      expected_by_dim[2] = TopologyInfo::faces( type );
    expected_by_dim[3] = elem ? 1 : 0;
    seen[2].resize( TopologyInfo::faces(type), false );
    seen[3].resize( 1, false );
  }
  
    // total number of samples
  const unsigned total = expected_by_dim[0] + expected_by_dim[1] 
                       + expected_by_dim[2] + expected_by_dim[3];
  
    // Check total sample count
  SamplePoints sp( corners, edges, faces, elem );
  CPPUNIT_ASSERT_EQUAL( total, sp.num_sample_points( type ) );
  

  for (unsigned i = 0; i < total; ++i)
  {
      // Check valid output from location_from_sample_number
    unsigned dim, num;
    sp.location_from_sample_number( type, i, dim, num );
    CPPUNIT_ASSERT( dim <= TopologyInfo::dimension( type ) );
    CPPUNIT_ASSERT( num < seen[dim].size() );
    
      // Check that sample_number_from_location is inverse of
      // location_from_sample_number
    unsigned sample = sp.sample_number_from_location( type, dim, num );
    CPPUNIT_ASSERT_EQUAL( i, sample );
    
      // Check that we don't get the same location more than once
    CPPUNIT_ASSERT_EQUAL( false, (bool)(seen[dim][num]) );
    seen[dim][num] = true;
    ++seen_by_dim[dim];
  }
  
    // Check that we got the correct sample points
  for (unsigned j = 0; j < 4; ++j)
  {
    CPPUNIT_ASSERT_EQUAL( expected_by_dim[j], seen_by_dim[j] );
  }
}
  
void SamplePointsTest::test_corners ()
{
  test_topo_common( true, false, false, false, TRIANGLE );
  test_topo_common( true, false, false, false, QUADRILATERAL );
  test_topo_common( true, false, false, false, TETRAHEDRON );
  test_topo_common( true, false, false, false, HEXAHEDRON );
}

void SamplePointsTest::test_edges ()
{
  test_topo_common( false, true, false, false, TRIANGLE );
  test_topo_common( false, true, false, false, QUADRILATERAL );
  test_topo_common( false, true, false, false, TETRAHEDRON );
  test_topo_common( false, true, false, false, HEXAHEDRON );
}

void SamplePointsTest::test_corners_and_edges ()
{
  test_topo_common( true, true, false, false, TRIANGLE );
  test_topo_common( true, true, false, false, QUADRILATERAL );
  test_topo_common( true, true, false, false, TETRAHEDRON );
  test_topo_common( true, true, false, false, HEXAHEDRON );
}

void SamplePointsTest::test_edges_and_faces ()
{
  test_topo_common( false, true, true, false, TRIANGLE );
  test_topo_common( false, true, true, false, QUADRILATERAL );
  test_topo_common( false, true, true, false, TETRAHEDRON );
  test_topo_common( false, true, true, false, HEXAHEDRON );
}

void SamplePointsTest::test_elem ()
{
  test_topo_common( false, false, false, true, PRISM );
  test_topo_common( false, false, false, true, PYRAMID );
  test_topo_common( false, false, false, true, TETRAHEDRON );
  test_topo_common( false, false, false, true, HEXAHEDRON );
}

void SamplePointsTest::test_corners_and_elem ()
{
  test_topo_common( true, false, false, true, PRISM );
  test_topo_common( true, false, false, true, PYRAMID );
  test_topo_common( true, false, false, true, TETRAHEDRON );
  test_topo_common( true, false, false, true, HEXAHEDRON );
}
