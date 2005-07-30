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
//  LAST-MOD:  9-Apr-03 at 10:05:06 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqVertexTest.cpp

Unit testing of various functions in the MsqVertex class. 

 */
// DESCRIP-END.
//


#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "PatchDataInstances.hpp"

#include <math.h>
#include <iostream>

#include "cppunit/extensions/HelperMacros.h"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class MsqVertexTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqVertexTest);
  CPPUNIT_TEST (test_flags);
  CPPUNIT_TEST (test_compute_weighted_jacobian_ideal_tri);
  CPPUNIT_TEST (test_compute_weighted_jacobian_ideal_quad);
  CPPUNIT_TEST (test_compute_weighted_jacobian_ideal_tet);
  CPPUNIT_TEST (test_compute_weighted_jacobian_ideal_hex);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData one_hex_patch;
  PatchData one_tet_patch;
  PatchData one_qua_patch;
  PatchData one_tri_patch;
  Vector3D e1, e2, e3;
  double tolEps;
  
public:
  void setUp()
  {
    tolEps=1.e-12;
      // set up the unit vectors
    e1.set(1,0,0);
    e2.set(0,1,0);
    e3.set(0,0,1);
    
    MsqPrintError err(cout);
    
    
    double hcoords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         2.0, 2.0, 1.0,
                         1.0, 2.0, 1.0,
                         1.0, 1.0, 2.0,
                         2.0, 1.0, 2.0,
                         2.0, 2.0, 2.0,
                         1.0, 2.0, 2.0 };
    size_t hconn[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    one_hex_patch.fill( 8, hcoords, 1, HEXAHEDRON, hconn, 0, err );
    

    double tcoords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         1.5, 1+sqrt(3.0)/2.0, 1.0,
                         1.5, 1+sqrt(3.0)/6.0, 1+sqrt(2.0)/sqrt(3.0) };
    size_t tconn[] = { 0, 1, 2, 3 };
    one_tet_patch.fill( 4, tcoords, 1, TETRAHEDRON, tconn, 0, err );

    double qcoords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         2.0, 2.0, 1.0,
                         1.0, 2.0, 1.0 };
    size_t qconn[] = { 0, 1, 2, 3 };
    one_qua_patch.fill( 4, qcoords, 1, QUADRILATERAL, qconn, 0, err );
    
    double rcoords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         1.5, 1+sqrt(3.0)/2.0, 1.0 };
    size_t rconn[] = { 0, 1, 2 };
    one_tri_patch.fill( 3, rcoords, 1, TRIANGLE, rconn, 0, err );
  }

  void tearDown()
  {

  }
  
public:
  MsqVertexTest()
    {}
  
  void test_hex_vertices()
  {
    MsqPrintError err(cout);
    // prints out the vertices.
    MsqVertex* ideal_vertices = one_hex_patch.get_vertex_array(err); CPPUNIT_ASSERT(!err);
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


  void test_compute_weighted_jacobian_ideal_hex()
  {
    MsqPrintError err(cout);
    //determinate of jacobs
    double deter;
    MsqMeshEntity* hex = one_hex_patch.get_element_array(err); CPPUNIT_ASSERT(!err);
    // get the ideal tets sample points 
    std::vector<Vector3D> sample_points;
    hex->get_sample_points(sample_points, err); CPPUNIT_ASSERT(!err);
      // and get the jacobian vectors (should be the identity matrix).
      //for the first sample point (0,0,0)
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    hex->compute_weighted_jacobian(one_hex_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); CPPUNIT_ASSERT(!err);
      //Make sure we have the identity
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[2]-e3).length()<tolEps);
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      hex->compute_weighted_jacobian(one_hex_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); CPPUNIT_ASSERT(!err);
        //make sure that our jacobian has determinat one
      deter=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
        //<<"\n\nInside test tet c_weighted_j deter = "<< deter <<"\n";
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
      //make sure that we had eight sample points.
    CPPUNIT_ASSERT(counter==8);
  }

       
  void test_compute_weighted_jacobian_ideal_tet()
  {
    MsqPrintError err(cout);
      //determinate of jacobs
    double deter;
    MsqMeshEntity* tet = one_tet_patch.get_element_array(err); CPPUNIT_ASSERT(!err);
      // get the ideal tets sample points 
    std::vector<Vector3D> sample_points;
    tet->get_sample_points(sample_points, err); CPPUNIT_ASSERT(!err);
    // and get the jacobian vectors (should be the identity matrix).
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    tet->compute_weighted_jacobian(one_tet_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); CPPUNIT_ASSERT(!err);

    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[2]-e3).length()<tolEps);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      tet->compute_weighted_jacobian(one_tet_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); CPPUNIT_ASSERT(!err);
      deter=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
        //<<"\n\nInside test tet c_weighted_j deter = "<< deter <<"\n";
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==4);
    
  }
void test_compute_weighted_jacobian_ideal_quad()
  {
    MsqPrintError err(cout);
      //determinate of jacobs
    double deter;
    MsqMeshEntity* quad = one_qua_patch.get_element_array(err); CPPUNIT_ASSERT(!err);
    // get the ideal quads sample points 
    std::vector<Vector3D> sample_points;
    quad->get_sample_points(sample_points, err); CPPUNIT_ASSERT(!err);
    // and get the jacobian vectors (should be the identity matrix).
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    quad->compute_weighted_jacobian(one_qua_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      quad->compute_weighted_jacobian(one_qua_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); CPPUNIT_ASSERT(!err);
        //deter is not the determinant in quad case
      deter=(jacobian_vectors[0]*jacobian_vectors[1]).length();
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==4);
    
  }
void test_compute_weighted_jacobian_ideal_tri()
  {
    MsqPrintError err(cout);
      //determinate of jacobs
    double deter;
    MsqMeshEntity* tri = one_tri_patch.get_element_array(err); CPPUNIT_ASSERT(!err);
    // get the ideal tris sample points 
    std::vector<Vector3D> sample_points;
    tri->get_sample_points(sample_points, err); CPPUNIT_ASSERT(!err);
    // and get the jacobian vectors (should be the identity matrix).
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    tri->compute_weighted_jacobian(one_tri_patch, *sp, jacobian_vectors,
                                   num_jac_vec , err); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      tri->compute_weighted_jacobian(one_tri_patch, *sp, jacobian_vectors,
                                     num_jac_vec , err); CPPUNIT_ASSERT(!err);
      deter=(jacobian_vectors[0]*jacobian_vectors[1]).length();
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==3);
  }

  void test_flags()
  {
     MsqVertex vtx(1,2,3);
     CPPUNIT_ASSERT( vtx.is_flag_set(MsqVertex::MSQ_HARD_FIXED)==false );
     vtx.set_hard_fixed_flag();
     CPPUNIT_ASSERT( vtx.is_flag_set(MsqVertex::MSQ_SOFT_FIXED)==false );
     CPPUNIT_ASSERT( vtx.is_flag_set(MsqVertex::MSQ_HARD_FIXED)==true );
     CPPUNIT_ASSERT( ( vtx.is_flag_set(MsqVertex::MSQ_SOFT_FIXED) ||
                       vtx.is_flag_set(MsqVertex::MSQ_HARD_FIXED) ) == true);
     CPPUNIT_ASSERT( ( vtx.is_flag_set(MsqVertex::MSQ_SOFT_FIXED) &&
                       vtx.is_flag_set(MsqVertex::MSQ_HARD_FIXED) ) == false);
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqVertexTest, "MsqVertexTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqVertexTest, "Unit");
