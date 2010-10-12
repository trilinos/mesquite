/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MeshUtilTest.cpp
 *  \brief Test MeshUtil class
 *  \author Jason Kraftcheck 
 */

#include "ArrayMesh.hpp"
#include "MeshUtil.hpp"
#include "UnitUtil.hpp"

using namespace Mesquite;

class MeshUtilTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(MeshUtilTest);
  CPPUNIT_TEST(test_edge_length_distribution_types);
  CPPUNIT_TEST(test_edge_length_distribution_unique);
  CPPUNIT_TEST(test_edge_length_distribution_empty);
  CPPUNIT_TEST_SUITE_END();

public:

    // test that edge_length_distribution function works
    // for all supported element types
  void test_edge_length_distribution_types();
  
    // test that edge_length_distribution function works
    // counts each implicit edge only once
  void test_edge_length_distribution_unique();
  
    // test that edge_length_distribution function 
    // handles the case of an empty mesh
  void test_edge_length_distribution_empty();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshUtilTest, "MeshUtilTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshUtilTest, "Unit");

const double EPSILON = 1e-8;

void MeshUtilTest::test_edge_length_distribution_types()
{
    // create a mesh that contains one of each element type,
    // completely disconnected from any other element
  double coords[] = { // triangle : 0 - 2
                      0, 0, 0,
                      2, 0, 0,
                      1, 1, 0,
                      // quad : 3 - 6
                      3, 0, 0,
                      3, 1, 0,
                      3, 2, 2,
                      3, 0, 1,
                      // tet : 7 - 10
                      0, 0, 0,
                      1, 0, 0,
                      0, 1, 0,
                      0, 0, 1,
                      // pyramid : 11 - 15
                      0, 0, 0,
                      0, -2, 0,
                      -2, -2, 0,
                      -1,  0, 0,
                      0, 0, -3,
                      // wedge : 16 - 21
                      0, 0, -2,
                      2, 0, -2,
                      1, 1, -2,
                      0, 0, -1.5,
                      2, 0, -1.5,
                      1, 1, -1,5,
                      // hex : 22 - 29
                      4, 0, 0,
                      4, 1, 0,
                      4, 2, 2,
                      4, 0, 1,
                      5, 0, 0,
                      5, 1, 0,
                      5, 2, 2,
                      5, 0, 1 };
  const int tri_offset = 0;
  const int quad_offset = tri_offset + 3;
  const int tet_offset = quad_offset + 4;
  const int pyr_offset = tet_offset + 4;
  const int wedge_offset = pyr_offset + 5;
  const int hex_offset = wedge_offset + 6;
  const int num_vtx = hex_offset + 8;
  const int fixed[num_vtx] = { 0 };
  EntityTopology types[] = { TRIANGLE, 
                             QUADRILATERAL,
                             TETRAHEDRON,
                             PYRAMID,
                             PRISM,
                             HEXAHEDRON };
  const unsigned long conn[] = { 0, 1, 2,     // triangle
                                 3, 4, 5, 6,  // quadrilateral
                                 7, 8, 9, 10, // tetrahedron
                                 11, 12, 13, 14, 15,            // pyramid
                                 16, 17, 18, 19, 20, 21,        // prism,
                                 22, 23, 24, 25, 26, 27, 28, 29 // hex 
                                };
  ArrayMesh mesh( 3, num_vtx, coords, fixed, 6, types, conn );


    // calculate expected value
  const int edges[][2] = { // triangle
                           { 0, 1 }, 
                           { 1, 2 },
                           { 2, 0 },
                           // quad
                           { 3, 4 },
                           { 4, 5 },
                           { 5, 6 },
                           { 6, 3 },
                           // tet
                           { 7, 8 },
                           { 8, 9 },
                           { 9, 7 },
                           { 7, 10 },
                           { 8, 10 },
                           { 9, 10 },
                           // pyramid
                           { 11, 12 },
                           { 12, 13 },
                           { 13, 14 },
                           { 14, 11 },
                           { 11, 15 },
                           { 12, 15 },
                           { 13, 15 },
                           { 14, 15 },
                           // prism
                           { 16, 17 },
                           { 17, 18 },
                           { 18, 16 },
                           { 16, 19 },
                           { 17, 20 },
                           { 18, 21 },
                           { 19, 20 },
                           { 20, 21 },
                           { 21, 19 },
                           // hex
                           { 22, 23 },
                           { 23, 24 },
                           { 24, 25 },
                           { 25, 22 },
                           { 22, 26 },
                           { 23, 27 },
                           { 24, 28 },
                           { 25, 29 },
                           { 26, 27 },
                           { 27, 28 },
                           { 28, 29 },
                           { 29, 26 } };
  const int num_edges = sizeof(edges)/sizeof(edges[0]);

  double exp_min = HUGE_VAL, exp_max = -HUGE_VAL;
  double exp_avg = 0, exp_rms = 0;
  for (int i = 0; i < num_edges; ++i) {
    int v1 = edges[i][0];
    int v2 = edges[i][1];
    Vector3D c0(coords + 3*v1);
    Vector3D c1(coords + 3*v2);
    Vector3D diff = c1 - c0;
    double len = sqrt(diff%diff);
    if (len < exp_min)
      exp_min = len;
    if (len > exp_max)
      exp_max = len;
    exp_avg += len;
    exp_rms += diff%diff;
  }
  
  exp_avg /= num_edges;
  exp_rms /= num_edges;
  double exp_std_dev = sqrt( exp_rms - exp_avg*exp_avg );
  exp_rms = sqrt(exp_rms);
  
  
  double act_min = -1, act_max = -1, act_avg = -1, act_rms = -1, act_std_dev = -1;
  MsqError err;
  MeshUtil tool(&mesh);
  tool.edge_length_distribution( act_min, act_avg, act_rms, act_max, act_std_dev, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_min, act_min, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_max, act_max, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_avg, act_avg, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_rms, act_rms, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_std_dev, act_std_dev, EPSILON );
}

void MeshUtilTest::test_edge_length_distribution_unique()
{
  // define two-tet mesh where tets share a face
  double coords[] = {  0, -5, 0,
                       0,  5, 0,
                       1,  0, 0,
                       0,  0, 0,
                       0,  0, 1 };
  const unsigned long conn[] = { 4, 3, 2, 0,
                                 2, 3, 4, 1 };
  int fixed[5] = {0};
  ArrayMesh mesh( 3, 5, coords, fixed, 2, TETRAHEDRON, conn );
  
  const int edges[][2] = { { 0, 2 },
                           { 0, 3 },
                           { 0, 4 },
                           { 1, 2 },
                           { 1, 3 },
                           { 1, 4 },
                           { 2, 3 },
                           { 3, 4 },
                           { 4, 2 } };
  const int num_edges = sizeof(edges)/sizeof(edges[0]);

  double exp_min = HUGE_VAL, exp_max = -HUGE_VAL;
  double exp_avg = 0, exp_rms = 0;
  for (int i = 0; i < num_edges; ++i) {
    int v1 = edges[i][0];
    int v2 = edges[i][1];
    Vector3D c0(coords + 3*v1);
    Vector3D c1(coords + 3*v2);
    Vector3D diff = c1 - c0;
    double len = sqrt(diff%diff);
    if (len < exp_min)
      exp_min = len;
    if (len > exp_max)
      exp_max = len;
    exp_avg += len;
    exp_rms += diff%diff;
  }
  
  exp_avg /= num_edges;
  exp_rms /= num_edges;
  double exp_std_dev = sqrt( exp_rms - exp_avg*exp_avg );
  exp_rms = sqrt(exp_rms);
  
  
  double act_min = -1, act_max = -1, act_avg = -1, act_rms = -1, act_std_dev = -1;
  MsqError err;
  MeshUtil tool(&mesh);
  tool.edge_length_distribution( act_min, act_avg, act_rms, act_max, act_std_dev, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_min, act_min, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_max, act_max, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_avg, act_avg, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_rms, act_rms, EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_std_dev, act_std_dev, EPSILON );
}

void MeshUtilTest::test_edge_length_distribution_empty()
{
  ArrayMesh empty( 3, 0, 0, 0, 0, 0, 0 );
  MeshUtil tool(&empty);
  double act_min = -1, act_max = -1, act_avg = -1, act_rms = -1, act_std_dev = -1;
  MsqError err;
  tool.edge_length_distribution( act_min, act_avg, act_rms, act_max, act_std_dev, err );
  CPPUNIT_ASSERT_EQUAL( MsqError::INVALID_MESH, err.error_code() );
}

