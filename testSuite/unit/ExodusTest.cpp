// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//

#include "Mesquite.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "MeshImpl.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

#include <iostream>

class ExodusTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(ExodusTest);
   CPPUNIT_TEST (test_elements);
   CPPUNIT_TEST_SUITE_END();

private:
  Mesquite::MeshImpl *mMesh;
  
public:
   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
#ifdef MSQ_USING_EXODUS
    Mesquite::MsqError err;
    
      // Read a Exodus Mesh file -- 10 triangles, 2 free vertices
    mMesh = new Mesquite::MeshImpl;
    mMesh->read_exodus("../../meshFiles/3D/CUBIT/test.g", err);
#endif
  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {
  }
  
public:
  ExodusTest()
    {}
  
  void test_elements()
  {
#ifdef MSQ_USING_EXODUS
    Mesquite::MsqError err;
    
      // Add mesh to a MeshSet.
    Mesquite::MeshSet mesh_set;
    mesh_set.add_mesh(mMesh, err); MSQ_CHKERR(err);

      // Get the number of vertices
    std::cout << "Number of vertices: "
              << mMesh->get_total_vertex_count(err) << std::endl;
      // Get the number of elements
    std::cout << "Number of elements: "
              << mMesh->get_total_element_count(err) << std::endl;
#endif
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ExodusTest, "ExodusTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ExodusTest, "Unit");
