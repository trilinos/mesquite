// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 22-Nov-02 at 08:05:56
//  LAST-MOD:  3-Dec-02 at 15:24:27 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MeshSetTest.cpp

Unit testing of various functions in the MeshSet class. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "MeshImpl.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

#include <list>
#include <iterator>

using namespace Mesquite;

class MeshSetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(MeshSetTest);
  CPPUNIT_WORK_IN_PROGRESS (test_add_multiple_meshes);
  CPPUNIT_TEST_SUITE_END();
  
private:
  
  Mesquite::MeshImpl *tri8;
  Mesquite::MeshImpl *quad4;
  
public:
    /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
    MsqError err;
    
      // Read a vtk file -- square meshed by 8 triangles
    tri8 = new Mesquite::MeshImpl;
    tri8->read_vtk("../../meshFiles/2D/VTK/square_tri_2.vtk", err);
    MSQ_CHKERR(err);
    
      // Read a vtk file
      // -- square meshed by 4 quads, adjacent to the previous mesh
    quad4 = new Mesquite::MeshImpl;
    quad4->read_vtk("../../meshFiles/2D/VTK/four_more_quads.vtk", err);
    MSQ_CHKERR(err);
  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {
  }
  
public:
  MeshSetTest()
   {}

   /*this function test the MeshSet concept of concatenating several meshes.*/
   void test_add_multiple_meshes()
   {
     MsqError err;
     
       /* Adds 2 adjacent meshes to the MeshSet. */
     MeshSet mesh_set;
     mesh_set.add_mesh(tri8, err); MSQ_CHKERR(err);
     mesh_set.add_mesh(quad4, err); MSQ_CHKERR(err);
     
       /* Retrieves a global patch */
     PatchData pd;
     PatchDataParameters pd_params;
     pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err, 0, 0);
     mesh_set.get_next_patch(pd, pd_params, err); MSQ_CHKERR(err);
     
     int num_vtx = pd.num_vertices();
     CPPUNIT_ASSERT( num_vtx == 18 );
   }
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshSetTest, "MeshSetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshSetTest, "Unit");
