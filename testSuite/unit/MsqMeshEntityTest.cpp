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
//  LAST-MOD: 13-Nov-02 at 16:48:45 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqMeshEntityTest.cpp

Unit testing of various functions in the MsqMeshEntity class. 

 */
// DESCRIP-END.
//


#include "MsqMeshEntity.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include <math.h>
#include <iostream>
#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;
using std::cout;

class MsqMeshEntityTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqMeshEntityTest);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_hex);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_tet);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_quad);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_tri);
  CPPUNIT_TEST (test_hex_vertices);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData one_hex_patch;
  PatchData one_tet_patch;
  PatchData one_qua_patch;
  PatchData one_tri_patch;
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
    
    coords[0] = 1.0; coords[1] = 1.0; coords[2] = 1.0;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2; coords[1] = 1; coords[2] = 1;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2.; coords[1] = 2.; coords[2] = 1;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 1.; coords[1] = 2.; coords[2] = 1;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 1.; coords[1] = 1.; coords[2] = 2;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2.; coords[1] = 1.; coords[2] = 2;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2.; coords[1] = 2.; coords[2] = 2;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 1.; coords[1] = 2.; coords[2] = 2;
    index = one_hex_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);
    
    // patch has only one element: an ideal hex.
    int indices[8];
    indices[0] = 0; indices[1] = 1; indices[2] = 2; indices[3] = 3;
    indices[4] = 4; indices[5] = 5; indices[6] = 6; indices[7] = 7;
    one_hex_patch.add_element(NULL, NULL, indices, HEXAHEDRON, err);
    MSQ_CHKERR(err);

      //**********************FILL TET*************************
      // creates empty Patch
    one_tet_patch.reserve_vertex_capacity (4, err); MSQ_CHKERR(err);
    one_tet_patch.reserve_element_capacity (1, err); MSQ_CHKERR(err);

    
      // Fills up with vertices for ideal tet
    
    coords[0] = 1; coords[1] = 1; coords[2] = 1;
    index = one_tet_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2; coords[1] = 1; coords[2] = 1;
    index = one_tet_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/2.0; coords[2] = 1;
    index = one_tet_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/6.0 ;
    coords[2] = 1+sqrt(2.0)/sqrt(3.0);
    index = one_tet_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);
    
      // patch has only one element: an ideal tet.
    int indices_tet[4];
    indices_tet[0] = 0; indices_tet[1] = 1; indices_tet[2] = 2;
    indices_tet[3] = 3;
    one_tet_patch.add_element(NULL, NULL, indices_tet, TETRAHEDRON, err);
    MSQ_CHKERR(err);
    
      //**********************FILL QUAD*************************
      // creates empty Patch
    one_qua_patch.reserve_vertex_capacity (4, err); MSQ_CHKERR(err);
    one_qua_patch.reserve_element_capacity (1, err); MSQ_CHKERR(err);

    
      // Fills up with vertices for ideal quad
    
    coords[0] = 1; coords[1] = 1; coords[2] = 1;
    index = one_qua_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2; coords[1] = 1; coords[2] = 1;
    index = one_qua_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2; coords[1] = 2; coords[2] = 1;
    index = one_qua_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 1; coords[1] = 2 ; coords[2] = 1;
    index = one_qua_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);
    
      // patch has only one element: an ideal quad.
    int indices_qua[4];
    indices_qua[0] = 0; indices_qua[1] = 1; indices_qua[2] = 2;
    indices_qua[3] = 3;
    one_qua_patch.add_element(NULL, NULL, indices_qua, QUADRILATERAL, err);
    MSQ_CHKERR(err);

          //**********************FILL tri*************************
      // creates empty Patch
    one_tri_patch.reserve_vertex_capacity (3, err); MSQ_CHKERR(err);
    one_tri_patch.reserve_element_capacity (1, err); MSQ_CHKERR(err);

    
      // Fills up with vertices for ideal tri
    
    coords[0] = 1; coords[1] = 1; coords[2] = 1;
    index = one_tri_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 2; coords[1] = 1; coords[2] = 1;
    index = one_tri_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);

    coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/2.0; coords[2] = 1;
    index = one_tri_patch.add_vertex(NULL, NULL, coords, false, err);
    MSQ_CHKERR(err);
    
      // patch has only one element: an ideal tri
    int indices_tri[3];
    indices_tri[0] = 0; indices_tri[1] = 1; indices_tri[2] = 2;
    one_tri_patch.add_element(NULL, NULL, indices_tri, TRIANGLE, err);
    MSQ_CHKERR(err);
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

    vtx.set(1,1,1);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[0]);
    
    vtx.set(2,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[6]);
    
    vtx.set(1,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[7]);
  }


  void test_compute_weigted_jacobian_ideal_hex()
  {
    MsqError err;
    //determinate of jacobs
    double deter;
    MsqMeshEntity* hex = one_hex_patch.get_element_array(err); MSQ_CHKERR(err);
    // get the ideal tets sample points 
    std::vector<Vector3D> sample_points;
    hex->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
      // and get the jacobian vectors (should be the identity matrix).
      //for the first sample point (0,0,0)
    int num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    hex->compute_weighted_jacobian(one_hex_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);
      //Make sure we have the identity
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<MSQ_MIN);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<MSQ_MIN);
    CPPUNIT_ASSERT((jacobian_vectors[2]-e3).length()<MSQ_MIN);
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      hex->compute_weighted_jacobian(one_hex_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
        //make sure that our jacobian has determinat one
      deter=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
        //<<"\n\nInside test tet c_weighted_j deter = "<< deter <<"\n";
      CPPUNIT_ASSERT(fabs(deter-1.0)<MSQ_MIN);
      ++sp;
      ++counter;
    }
      //make sure that we had eight sample points.
    CPPUNIT_ASSERT(counter==8);
  }

       
  void test_compute_weigted_jacobian_ideal_tet()
  {
    MsqError err;
      //determinate of jacobs
    double deter;
    MsqMeshEntity* tet = one_tet_patch.get_element_array(err); MSQ_CHKERR(err);
      // get the ideal tets sample points 
    std::vector<Vector3D> sample_points;
    tet->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
    // and get the jacobian vectors (should be the identity matrix).
    int num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    tet->compute_weighted_jacobian(one_tet_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);

    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<MSQ_MIN);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<MSQ_MIN);
    CPPUNIT_ASSERT((jacobian_vectors[2]-e3).length()<MSQ_MIN);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      tet->compute_weighted_jacobian(one_tet_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
      deter=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
        //<<"\n\nInside test tet c_weighted_j deter = "<< deter <<"\n";
      CPPUNIT_ASSERT(fabs(deter-1.0)<MSQ_MIN);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==4);
    
  }
void test_compute_weigted_jacobian_ideal_quad()
  {
    MsqError err;
      //determinate of jacobs
    double deter;
    MsqMeshEntity* quad = one_qua_patch.get_element_array(err); MSQ_CHKERR(err);
    // get the ideal quads sample points 
    std::vector<Vector3D> sample_points;
    quad->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
    // and get the jacobian vectors (should be the identity matrix).
    int num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    quad->compute_weighted_jacobian(one_qua_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<MSQ_MIN);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<MSQ_MIN);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      quad->compute_weighted_jacobian(one_qua_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
        //deter is not the determinant in quad case
      deter=(jacobian_vectors[0]*jacobian_vectors[1]).length();
      CPPUNIT_ASSERT(fabs(deter-1.0)<MSQ_MIN);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==4);
    
  }
void test_compute_weigted_jacobian_ideal_tri()
  {
    MsqError err;
      //determinate of jacobs
    double deter;
    MsqMeshEntity* tri = one_tri_patch.get_element_array(err); MSQ_CHKERR(err);
    // get the ideal tris sample points 
    std::vector<Vector3D> sample_points;
    tri->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
    // and get the jacobian vectors (should be the identity matrix).
    int num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    tri->compute_weighted_jacobian(one_tri_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<MSQ_MIN);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<MSQ_MIN);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      tri->compute_weighted_jacobian(one_tri_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
      deter=(jacobian_vectors[0]*jacobian_vectors[1]).length();
      CPPUNIT_ASSERT(fabs(deter-1.0)<MSQ_MIN);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==3);
    
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "MsqMeshEntityTest");
