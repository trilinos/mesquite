// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 14-Jan-02 at 08:05:56
//  LAST-MOD: 20-Jan-03 at 12:59:29 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file AomdVtkTest.cpp

Unit testing of the uploading of VTK format into AOMD.

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

class AomdVtkTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(AomdVtkTest);
   CPPUNIT_TEST (test_elements);
   CPPUNIT_TEST_SUITE_END();
   
private:
   
   TSTT::Mesh_Handle tri10;

public:
   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
     MsqError err;
     
     char file_name[128];
     TSTT::MeshError tstt_err;
     
     /* Reads a TSTT Mesh file -- 10 triangles, 2 free vertices */
     TSTT::Mesh_Create(&tri10, &tstt_err); assert(!tstt_err);
     strcpy(file_name, "../../meshFiles/2D/VTK/tri_10.vtk");
     std::cout << file_name << std::endl; // dbg
     TSTT::Mesh_Load(tri10, file_name, &tstt_err);
     std::cout << tstt_err << std::endl; // dbg
     assert(!tstt_err);
     
  }

   /* Automatically called by CppUnit after each test function. */
  void tearDown()
  {
  }
  
public:
  AomdVtkTest()
   {}

   /*  */
   void test_elements()
   {
      MsqError err;
     /* Adds TSTT mesh to a MeshSet. */
      MeshSet mesh_set;
      mesh_set.add_mesh(tri10, err); MSQ_CHKERR(err);

      /* Retrieves a global patch */
      PatchData pd;
      PatchDataParameters pd_params;
      pd_params.set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err, 1, 0);
      pd_params.add_culling_method(PatchData::NO_BOUNDARY_VTX);
      
      mesh_set.get_next_patch(pd, pd_params, err); MSQ_CHKERR(err);

      int free_vtx = pd.num_free_vertices(err); MSQ_CHKERR(err);
      std::cout << "nb of free vertices: " << free_vtx << std::endl;
      CPPUNIT_ASSERT( free_vtx == 1 );
      
      MsqMeshEntity* element_array =  pd.get_element_array(err); MSQ_CHKERR(err);
      int num_elements = pd.num_elements();
      CPPUNIT_ASSERT( num_elements == 6 );

      for (int i=0; i<num_elements; ++i) {
         std::cout << element_array[i];
      }
      
      MsqVertex* vtx_array = pd.get_vertex_array(err); MSQ_CHKERR(err);
      int num_vertices = pd.num_vertices();
      CPPUNIT_ASSERT( num_vertices == 7 );

      for (int i=0; i<num_vertices; ++i) {
         std::cout << vtx_array[i];
      }

      
      mesh_set.get_next_patch(pd, pd_params, err); MSQ_CHKERR(err);

      element_array =  pd.get_element_array(err); MSQ_CHKERR(err);
      num_elements = pd.num_elements();
      CPPUNIT_ASSERT( num_elements == 6 );

      for (int i=0; i<num_elements; ++i) {
         std::cout << element_array[i];
      }
      
      vtx_array = pd.get_vertex_array(err); MSQ_CHKERR(err);
      num_vertices = pd.num_vertices();
      CPPUNIT_ASSERT( num_vertices == 7 );

      for (int i=0; i<num_vertices; ++i) {
         std::cout << vtx_array[i];
      }
   }
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AomdVtkTest, "AomdVtkTest");

