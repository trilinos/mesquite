// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
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

#include "TSTT_Base.h"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

#include <list>
#include <iterator>

using namespace Mesquite;

class MeshSetTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(MeshSetTest);
   CPPUNIT_TEST (test_add_multiple_meshes);
   CPPUNIT_TEST_SUITE_END();
   
private:
   
   TSTT::Mesh_Handle tri8;
   TSTT::Mesh_Handle quad4;

public:
   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
     MsqError err;
     
     char file_name[128];
     TSTT::MeshError tstt_err;
     
     /* Reads a TSTT Mesh file -- square meshed by 8 triangles */
     TSTT::Mesh_Create(&tri8, &tstt_err); assert(!tstt_err);
     strcpy(file_name, "../../meshFiles/2D/VTK/square_tri_2.vtk");
     std::cout << file_name << std::endl; // dbg
     TSTT::Mesh_Load(tri8, file_name, &tstt_err);
     std::cout << tstt_err << std::endl; // dbg
     assert(!tstt_err);
     
     /* Reads a TSTT Mesh file -- square meshed by 4 quads, adjacent to the previous mesh */
     TSTT::Mesh_Create(&quad4, &tstt_err); assert(!tstt_err);
     strcpy(file_name, "../../meshFiles/2D/VTK/four_more_quads.vtk");
     std::cout << file_name << std::endl; // dbg
     TSTT::Mesh_Load(quad4, file_name, &tstt_err);
     std::cout << tstt_err << std::endl; // dbg
     assert(!tstt_err);
  }

   /* Automatically called by CppUnit after each test function. */
  void tearDown()
  {
  }
  
public:
  MeshSetTest()
   {}

   /* this function test the MeshSet concept of concatenating several TSTT meshes. */
   void test_add_multiple_meshes()
   {
      MsqError err;

      /* Adds 2 adjacent TSTT meshes to the MeshSet. */
      MeshSet mesh_set;
      mesh_set.add_mesh(tri8, err); MSQ_CHKERR(err);
      mesh_set.add_mesh(quad4, err); MSQ_CHKERR(err);

      /* Retrieves a global patch */
      PatchData pd;
      PatchDataParameters pd_params;
      pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err, 0, 0);
      mesh_set.get_next_patch(pd, pd_params, err); MSQ_CHKERR(err);

      int num_vtx = pd.num_vertices();
      CPPUNIT_ASSERT( num_vtx == 15 );

      
   }
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshSetTest, "MeshSetTest");

