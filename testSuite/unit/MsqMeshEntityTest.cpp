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
//  LAST-MOD:  5-Dec-02 at 16:55:01 by Thomas Leurent
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
#include "PatchDataInstances.hpp"
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
    create_one_hex_patch(one_hex_patch, err); MSQ_CHKERR(err);
    create_one_tet_patch(one_tet_patch, err); MSQ_CHKERR(err);
    create_one_tri_patch(one_tri_patch, err); MSQ_CHKERR(err);
    create_one_quad_patch(one_qua_patch, err); MSQ_CHKERR(err);
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
