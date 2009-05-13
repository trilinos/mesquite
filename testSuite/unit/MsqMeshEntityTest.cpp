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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 12-Nov-02 at 18:05:56
//  LAST-MOD:  7-May-03 at 13:28:47 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqMeshEntityTest.cpp

Unit testing of various functions in the MsqMeshEntity class.

\author Michael Brewer
\author Thomas Leurent

 */
// DESCRIP-END.
//


#include "MsqMeshEntity.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "PatchDataInstances.hpp"
#include <math.h>
#include <iostream>
#include "UnitUtil.hpp"

using namespace Mesquite;
using std::cout;
using std::endl;

class MsqMeshEntityTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqMeshEntityTest);
  CPPUNIT_TEST (test_hex_vertices);
  CPPUNIT_TEST (test_centroid_tri);
  CPPUNIT_TEST (test_centroid_quad);
  CPPUNIT_TEST (test_centroid_hex);
  CPPUNIT_TEST (test_unsigned_area);
  CPPUNIT_TEST (test_unsigned_area_poly);
  CPPUNIT_TEST (test_unsigned_area_tet);
  CPPUNIT_TEST (test_unsigned_area_pyr);
  CPPUNIT_TEST (test_unsigned_area_pri);
  CPPUNIT_TEST (test_unsigned_area_hex);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData oneHexPatch;
  PatchData oneTetPatch;
  PatchData oneQuadPatch;
  PatchData oneTriPatch;
  Vector3D e1, e2, e3;
  double tolEps;
public:
  void setUp()
  {
    tolEps=1.e-12;
    
    // sets up the unit vectors
    e1.set(1,0,0);
    e2.set(0,1,0);
    e3.set(0,0,1);

    MsqPrintError err(cout);

    // creates empty Patch
    create_one_hex_patch(oneHexPatch, err); CPPUNIT_ASSERT(!err);
    create_one_tet_patch(oneTetPatch, err); CPPUNIT_ASSERT(!err);
    create_one_tri_patch(oneTriPatch, err); CPPUNIT_ASSERT(!err);
    create_one_quad_patch(oneQuadPatch, err); CPPUNIT_ASSERT(!err);
  }

  void tearDown()
  {
    destroy_patch_with_domain( oneTriPatch );
    destroy_patch_with_domain( oneQuadPatch );
  }
  
public:
  MsqMeshEntityTest()
  {  }
  
  void test_hex_vertices()
  {
    MsqPrintError err(cout);
    // prints out the vertices.
    const MsqVertex* ideal_vertices = oneHexPatch.get_vertex_array(err); CPPUNIT_ASSERT(!err);
    size_t num_vtx = oneHexPatch.num_nodes();
    CPPUNIT_ASSERT_EQUAL(size_t(8), num_vtx);
    
    MsqVertex vtx;

    vtx.set(1,1,1);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[0]);
    
    vtx.set(2,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[6]);
    
    vtx.set(1,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[7]);
  }

  //! test the centroid of the first element in the Patch
  void test_centroid(PatchData &pd, Vector3D &correct)
  {
    MsqPrintError err(cout);
    double eps = 1e-6;
    Vector3D centroid;

    MsqMeshEntity* elem = pd.get_element_array(err); CPPUNIT_ASSERT(!err);
    elem->get_centroid(centroid, pd, err); CPPUNIT_ASSERT(!err);

//     cout << "centroid: "<< centroid <<endl; 
//     cout << "correct: "<< correct <<endl; 

    for (int i=0; i<3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL( centroid[i], correct[i], eps);
  }

  void test_centroid_tri()
  {
    Vector3D correct(1.5, 1+1/(2.0*sqrt(3.0)), 1.0);    
    test_centroid(oneTriPatch, correct);
  }

  void test_centroid_quad()
  {
    Vector3D correct(1.5, 1.5, 1.0);    
    test_centroid(oneQuadPatch, correct);
  }

  void test_centroid_hex()
  {
    Vector3D correct(1.5, 1.5, 1.5);    
    test_centroid(oneHexPatch, correct);
  }


  void test_unsigned_area()
     {
       MsqPrintError err(cout);
       MsqMeshEntity* tri = oneTriPatch.get_element_array(err);
       CPPUNIT_ASSERT(!err);
       CPPUNIT_ASSERT(fabs(tri->compute_unsigned_area(oneTriPatch,err)
                           -(sqrt(3.0)/4.0)) < tolEps);
       MsqMeshEntity* quad = oneQuadPatch.get_element_array(err);
       CPPUNIT_ASSERT(!err);
       CPPUNIT_ASSERT(fabs(quad->compute_unsigned_area(oneQuadPatch,err)
                           -1.0) < tolEps);
     }

  void test_unsigned_area_poly();
  void test_unsigned_area_tet();
  void test_unsigned_area_pyr();
  void test_unsigned_area_pri();
  void test_unsigned_area_hex();
  
  void test_unsigned_area_common( EntityTopology type,
                                  const double* coords,
                                  double expected );
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "MsqMeshEntityTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "Unit");
 
const size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
const bool fixed[] = { false, false, false, false, false, false, false, false };

void MsqMeshEntityTest::test_unsigned_area_poly()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                          0.5, 1.5, 0,
                            0, 1, 0 };
  size_t n_vtx = 5;
  EntityTopology type = POLYGON;
  MsqError err;
  
  PatchData pd;
  pd.fill( n_vtx, coords, 1, &type, &n_vtx, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  double a = pd.element_by_index(0).compute_unsigned_area( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.25, a, 1e-8 );
}

void MsqMeshEntityTest::test_unsigned_area_common( EntityTopology type,
                                                   const double* coords,
                                                   double expected )
{
  MsqError err;
  PatchData pd;

  pd.fill( TopologyInfo::corners( type ), coords, 1, type, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  double a = pd.element_by_index(0).compute_unsigned_area( pd, err );
  ASSERT_NO_ERROR(err);
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, a, 1e-8 );
}


void MsqMeshEntityTest::test_unsigned_area_tet()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            0, 1, 0,
                            0, 0, 1 };
  test_unsigned_area_common( TETRAHEDRON, coords, 1.0/6.0 );
}

void MsqMeshEntityTest::test_unsigned_area_pyr()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                            0, 1, 0,
                            0, 0, 1 };
  test_unsigned_area_common( PYRAMID, coords, 1.0/3.0 );

  const double pyr_coords[] = {-1, -1, -1,
                                1, -1, -1,
                                1,  1, -1,
                               -1,  1, -1,
                                0,  0,  0 };
  test_unsigned_area_common( PYRAMID, pyr_coords, 4.0/3.0 );
}

void MsqMeshEntityTest::test_unsigned_area_pri()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            0, 1, 0,
                            0, 0, 1,
                            1, 0, 1,
                            0, 1, 1 };
  test_unsigned_area_common( PRISM, coords, 0.5 );
  
  const double tet_coords[] = { 0, 0, 0,
                                1, 0, 0,
                                0, 1, 0,
                                0, 0, 1,
                                0, 0, 1,
                                0, 0, 1 };
  test_unsigned_area_common( PRISM, tet_coords, 1.0/6.0 );
}

void MsqMeshEntityTest::test_unsigned_area_hex()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                            0, 1, 0,
                            0, 0, 1,
                            1, 0, 1,
                            1, 1, 1,
                            0, 1, 1 };
  test_unsigned_area_common( HEXAHEDRON, coords, 1.0 );

  const double coords2[] = { 0, 0, 0,
                             2, 0, 0,
                             2, 2, 0,
                             0, 2, 0,
                             0, 0, 2,
                             2, 0, 2,
                             2, 2, 2,
                             0, 2, 2 };
  test_unsigned_area_common( HEXAHEDRON, coords2, 8.0 );
  
  const double pyr_coords[] = {-1,-1, 0,
                                1,-1, 0,
                                1, 1, 0,
                               -1, 1, 0,
                                0, 0, 1,
                                0, 0, 1,
                                0, 0, 1,
                                0, 0, 1 };
  test_unsigned_area_common( HEXAHEDRON, pyr_coords, 4.0/3.0 );
}


