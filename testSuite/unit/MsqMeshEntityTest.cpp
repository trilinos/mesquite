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
//  LAST-MOD: 13-Nov-02 at 15:44:35 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqMeshEntityTest.cpp

Unit testing of various functions in the MsqMeshEntity class. 

 */
// DESCRIP-END.
//


// TODO : get rid off ? 
#ifndef MsqMeshEntityTest_cpp
#define MsqMeshEntityTest_cpp

#include "MsqMeshEntity.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;

class MsqMeshEntityTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqMeshEntityTest);
  CPPUNIT_TEST (test_compute_weigthed_jacobian_ideal_hex);
  CPPUNIT_TEST (test_hex_vertices);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData one_hex_patch;
  Vector3D e1, e2, e3;

public:
  void setUp()
  {
    // sets up the unit vectors
    e1.set(1,0,0);
    e2.set(0,1,0);
    e3.set(0,0,1);

    MsqError err;

    // creates empty Patch
    one_hex_patch.reserve_vertex_capacity (8, err); MSQ_CHKERR(err);
    one_hex_patch.reserve_element_capacity (1, err); MSQ_CHKERR(err);

    
    // Fills up with vertices for ideal hexahedra.
    double coords[3];
    int index;
    
    coords[0] = 0; coords[1] = 0; coords[2] = 0;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    coords[0] = 1; coords[1] = 0; coords[2] = 0;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    coords[0] = 1.; coords[1] = 1.; coords[2] = 0;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    coords[0] = 0.; coords[1] = 1.; coords[2] = 0;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    coords[0] = 0.; coords[1] = 0.; coords[2] = 1;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    coords[0] = 1.; coords[1] = 0.; coords[2] = 1;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    coords[0] = 1.; coords[1] = 1.; coords[2] = 1;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    coords[0] = 0.; coords[1] = 1.; coords[2] = 1;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err); MSQ_CHKERR(err);

    
    // patch has only one element: an ideal tet.
    int indices[8];
    indices[0] = 0; indices[1] = 1; indices[2] = 2; indices[3] = 3;
    indices[4] = 4; indices[5] = 5; indices[6] = 6; indices[7] = 7;
    one_hex_patch.add_element(NULL, NULL, indices, HEXAHEDRON, err); MSQ_CHKERR(err);

  }

  void tearDown()
  {

  }
  
public:
  MsqMeshEntityTest()
    {}
  
  void test_hex_vertices()
  {
    MsqError err;
    // prints out the vertices.
    MsqVertex* ideal_vertices = one_hex_patch.get_vertex_array(err); MSQ_CHKERR(err);
    int num_vtx = one_hex_patch.num_vertices();
    CPPUNIT_ASSERT_EQUAL(8, num_vtx);
    
    MsqVertex vtx;

    vtx.set(0,0,0);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[0]);
    
    vtx.set(1,1,1);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[6]);
    
    vtx.set(0,1,1);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[7]);
  }


  void test_compute_weigthed_jacobian_ideal_hex()
  {
    MsqError err;
    MsqMeshEntity* hex = one_hex_patch.get_element_array(err); MSQ_CHKERR(err);
    // get the ideal tets sample points 
    std::vector<Vector3D> sample_points;
    hex->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
    // and get the jacobian vectors (should be the identity matrix).
    int num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    hex->compute_weighted_jacobian(one_hex_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);

    CPPUNIT_ASSERT(jacobian_vectors[0] == e1);
    CPPUNIT_ASSERT(jacobian_vectors[1] == e2);
    CPPUNIT_ASSERT(jacobian_vectors[2] == e3);
  }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "MsqMeshEntityTest");

#endif
