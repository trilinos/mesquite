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
//  LAST-MOD: 12-Nov-02 at 18:10:15 by Thomas Leurent
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
#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

class MsqMeshEntityTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqMeshEntityTest);
  CPPUNIT_TEST (test_default_constructor);
  CPPUNIT_TEST (test_double_constructor);
  CPPUNIT_TEST (test_copy_constructor);
  CPPUNIT_TEST_EXCEPTION (throw_exception, CppUnit::SignalException);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData one_hex_patch;
  Vector3D e1, e2, e3;

protected:
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
    coords[0] = 0; coords[1] = 0; coords[2] = 0;
    cout << one_hex_patch.add_vertex(coords, false, err); MSQ_CHKERR(err);
    coords[0] = 1; coords[1] = 0; coords[2] = 0;
    cout << one_hex_patch.add_vertex(coords, false, err); MSQ_CHKERR(err);
    coords[0] = 1.; coords[1] = 1.; coords[2] = 0;
    cout << one_hex_patch.add_vertex(coords, false, err); MSQ_CHKERR(err);
    coords[0] = 0.; coords[1] = 1.; coords[2] = 0;
    cout << one_hex_patch.add_vertex(coords, false, err) << endl; MSQ_CHKERR(err);
    coords[0] = 0.; coords[1] = 0.; coords[2] = 1;
    cout << one_hex_patch.add_vertex(coords, false, err) << endl; MSQ_CHKERR(err);
    coords[0] = 1.; coords[1] = 0.; coords[2] = 1;
    cout << one_hex_patch.add_vertex(coords, false, err) << endl; MSQ_CHKERR(err);
    coords[0] = 1.; coords[1] = 1.; coords[2] = 1;
    cout << one_hex_patch.add_vertex(coords, false, err) << endl; MSQ_CHKERR(err);
    coords[0] = 0.; coords[1] = 1.; coords[2] = 1;
    cout << one_hex_patch.add_vertex(coords, false, err) << endl; MSQ_CHKERR(err);
    // patch has only one element: an ideal tet.
    int indices[8];
    indices[0] = 0; indices[1] = 1; indices[2] = 2; indices[3] = 3;
    indices[4] = 4; indices[5] = 5; indices[6] = 6; indices[7] = 7;
    one_hex_patch.add_element(indices, HEXAHEDRON, me); MSQ_CHKERR(me);

    // prints out the vertices.
    MsqVertex* ideal_vertices = simple_patch.get_vertex_array(me); MSQ_CHKERR(me);
    int num_vtx = simple_patch.num_vertices();
    cout << endl << "ideal Hex vertices \n";
    for (int i=0; i<num_vtx; ++i)
      cout << ideal_vertices[i];
  }

  void tearDown()
  {

  }
  
public:
  MsqMeshEntityTest()
    {}
  
//   void test_default_constructor()
//     {
//     }

  void test_compute_weigthed_jacobian_ideal_hex()
    {
      MsqMeshEntity* hex = one_hex_patch.get_element_array(me); MSQ_CHKERR(me);
      // get the ideal tets sample points 
      std::vector<Vector3D> sample_points;
      hex->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, me); MSQ_CHKERR(err);
      // prints the sample points
      std::vector<Vector3D>::iterator sp;
      cout << endl << "Hexahedra sample points \n";
      for (sp=sample_points.begin(); sp!=sample_points.end(); ++sp)
        cout << *sp;
      // and get the jacobian vectors (should be the identity matrix).
      int num_jac_vec;
      sp=sample_points.begin();
      Vector3D jacobian_vectors[3];
      hex->compute_weighted_jacobian(simple_patch, &(*sp), jacobian_vectors,
                                     num_jac_vec , me); MSQ_CHKERR(me);
      cout << endl << "Jacobian matrix \n";
      for (int i=0; i< num_jac_vec; ++i)
        cout << jacobian_vectors[i];

      CPPUNIT_ASSERT(jacobian_vectors[0] == e1);
      CPPUNIT_ASSERT(jacobian_vectors[1] == e2);
      CPPUNIT_ASSERT(jacobian_vectors[2] == e3);
    }

//   void test_copy_constructor()
//     {
//     }

//   void test_equals()
//     {
//       CPPUNIT_ASSERT(v1 == v2);
//     }

//   void throw_exception()
//     {
//       Mesquite::Vector3D* v = NULL;
//       double d = v->x();
//       std::cout << d << std::endl;
//     }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "Misc");

#endif
