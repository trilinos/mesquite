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
//  LAST-MOD: 14-Jan-03 at 11:34:05 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PatchDataInstances.hpp

This header file contains some functions to instantiates particular PatchData Objects.
Those objects can be used in unit tests.
Patches must be allocated and dealocated by the caller. 

*/
// DESCRIP-END.
//

#ifndef PatchDataInstances_hpp
#define PatchDataInstances_hpp

#include "MsqVertex.hpp"
#include "PatchData.hpp"

#include <math.h>
#include <iostream>

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

namespace Mesquite {

   //! creates a patch containing one ideal tetrahedra
   inline void create_one_hex_patch(PatchData &one_hex_patch, MsqError &err) {
      
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
      size_t indices[8];
      indices[0] = 0; indices[1] = 1; indices[2] = 2; indices[3] = 3;
      indices[4] = 4; indices[5] = 5; indices[6] = 6; indices[7] = 7;
      one_hex_patch.add_element(NULL, NULL, indices, HEXAHEDRON, err);
      MSQ_CHKERR(err);
   }


   //! creates a Patch containing an ideal tetrahedra
   inline void create_one_tet_patch(PatchData &one_tet_patch, MsqError &err) {
      //**********************FILL TET*************************
      // creates empty Patch
      one_tet_patch.reserve_vertex_capacity (4, err); MSQ_CHKERR(err);
      one_tet_patch.reserve_element_capacity (1, err); MSQ_CHKERR(err);

      // Fills up with vertices for ideal tet
      double coords[3];
      int index;
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
      size_t indices_tet[4];
      indices_tet[0] = 0; indices_tet[1] = 1; indices_tet[2] = 2;
      indices_tet[3] = 3;
      one_tet_patch.add_element(NULL, NULL, indices_tet, TETRAHEDRON, err);
      MSQ_CHKERR(err);
   }
    

   //! creates a Patch containing an ideal quadrilateral
   inline void create_one_quad_patch(PatchData &one_qua_patch, MsqError &err) {
      //**********************FILL QUAD*************************
      // creates empty Patch
      one_qua_patch.reserve_vertex_capacity (4, err); MSQ_CHKERR(err);
      one_qua_patch.reserve_element_capacity (1, err); MSQ_CHKERR(err);
    
      // Fills up with vertices for ideal quad
      double coords[3];
      int index;
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
      size_t indices_qua[4];
      indices_qua[0] = 0; indices_qua[1] = 1; indices_qua[2] = 2;
      indices_qua[3] = 3;
      one_qua_patch.add_element(NULL, NULL, indices_qua, QUADRILATERAL, err);
      MSQ_CHKERR(err);
   }


   /*! \fn create_one_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing an ideal triangle
          /   \
         0-----1 
   */
   inline void create_one_tri_patch(PatchData &one_tri_patch, MsqError &err) {
      //**********************FILL tri*************************
      // creates empty Patch
      one_tri_patch.reserve_vertex_capacity (3, err); MSQ_CHKERR(err);
      one_tri_patch.reserve_element_capacity (1, err); MSQ_CHKERR(err);
    
      // Fills up with vertices for ideal tri
      double coords[3];
      int index;
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
      size_t indices_tri[3];
      indices_tri[0] = 0; indices_tri[1] = 1; indices_tri[2] = 2;
      one_tri_patch.add_element(NULL, NULL, indices_tri, TRIANGLE, err);
      MSQ_CHKERR(err);
   }

  
   /*! \fn create_two_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing two ideal triangles
          / 0 \
         0-----1
          \ 1 /
           \ /
            3
   */
   inline void create_two_tri_patch(PatchData &pd, MsqError &err) {
      //**********************FILL tri*************************
      // creates empty Patch
      pd.reserve_vertex_capacity (4, err); MSQ_CHKERR(err);
      pd.reserve_element_capacity (2, err); MSQ_CHKERR(err);
    
      // Fills up with vertices for ideal triangles
      double coords[3];
      int index;
      coords[0] = 1; coords[1] = 1; coords[2] = 1;
      index = pd.add_vertex(NULL, NULL, coords, false, err);
      MSQ_CHKERR(err);
      coords[0] = 2; coords[1] = 1; coords[2] = 1;
      index = pd.add_vertex(NULL, NULL, coords, false, err);
      MSQ_CHKERR(err);
      coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/2.0; coords[2] = 1;
      index = pd.add_vertex(NULL, NULL, coords, false, err);
      MSQ_CHKERR(err);
      coords[0] = 1.5; coords[1] = 1-sqrt(3.0)/2.0; coords[2] = 1;
      index = pd.add_vertex(NULL, NULL, coords, false, err);
      MSQ_CHKERR(err);
    
      // patch has only one element: an ideal tri
      size_t indices_tri[3];
      indices_tri[0] = 0; indices_tri[1] = 1; indices_tri[2] = 2;
      pd.add_element(NULL, NULL, indices_tri, TRIANGLE, err);
      MSQ_CHKERR(err);
      indices_tri[0] = 0; indices_tri[1] = 3; indices_tri[2] = 1;
      pd.add_element(NULL, NULL, indices_tri, TRIANGLE, err);
      MSQ_CHKERR(err);
   }
  

   /*! \fn create_four_quads_patch(PatchData &four_quads, MsqError &err)
     our 2D set up: 4 quads, center vertex outcentered by (0,-0.5)
      7____6____5
      |    |    |
      | 2  |  3 |
      8-_  |  _-4       vertex 1 is at (0,0)
      |  -_0_-  |       vertex 5 is at (2,2)
      | 0  |  1 |
      1----2----3
   */
   inline void create_four_quads_patch(PatchData &four_quads, MsqError &err) 
   {
      four_quads.reserve_vertex_capacity(9, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 1,.5, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 0, 0, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 1, 0, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 2, 0, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 2, 1, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 2, 2, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 1, 2, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 0, 2, 0, true, err); MSQ_CHKERR(err);
      four_quads.add_vertex(NULL, NULL, 0, 1, 0, true, err); MSQ_CHKERR(err);
      
      size_t ind[4];
      four_quads.reserve_element_capacity(4, err); MSQ_CHKERR(err);
      ind[0] = 1; ind[1]=2; ind[2]=0; ind[3]=8;
      four_quads.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 2; ind[1]=3; ind[2]=4; ind[3]=0;
      four_quads.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 8; ind[1]=0; ind[2]=6; ind[3]=7;
      four_quads.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 0; ind[1]=4; ind[2]=5; ind[3]=6;
      four_quads.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
   }
   

   /*! \fn create_six_quads_patch(PatchData &four_quads, MsqError &err)
     our 2D set up: 6 quads, 1 center vertex outcentered by (0,-0.5), the other centered
      7____6____5___11
      |    |    |    |
      | 2  |  3 | 5  |
      8-_  |  _-4---10       vertex 1 is at (0,0)
      |  -_0_-  |    |       vertex 11 is at (3,2)
      | 0  |  1 | 4  |
      1----2----3----9
   */
   inline void create_six_quads_patch(PatchData &pd, MsqError &err) 
   {
      pd.reserve_vertex_capacity(12, err); MSQ_CHKERR(err);
      pd.add_vertex(NULL, NULL, 1,.5, 0, true, err, MsqVertex::MSQ_NO_VTX_FLAG); 
      pd.add_vertex(NULL, NULL, 0, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 1, 0, true, err, MsqVertex::MSQ_NO_VTX_FLAG); 
      pd.add_vertex(NULL, NULL, 2, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 1, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 1, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      
      size_t ind[4];
      pd.reserve_element_capacity(6, err); MSQ_CHKERR(err);
      ind[0] = 1; ind[1]=2; ind[2]=0; ind[3]=8;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 2; ind[1]=3; ind[2]=4; ind[3]=0;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 8; ind[1]=0; ind[2]=6; ind[3]=7;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 0; ind[1]=4; ind[2]=5; ind[3]=6;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 3; ind[1]=9; ind[2]=10; ind[3]=4;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 4; ind[1]=10; ind[2]=11; ind[3]=5;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
   }


   /*! \fn create_twelve_hex_patch(PatchData &pd, MsqError &err)
     3D set up: 12 quads, one center vertex outcentered by (0,-0.5),
     the other centered. Vertex 1 is at (0,0,-1). Vertex 35 is at (3,2,1).
     
      7____6____5___11     19___18____17__23     31___30___29___35
      |    |    |    |      |    |    |    |      |    |    |    |
      | 2  |  3 | 5  |      |    |    |    |      | 8  |  9 | 11 |
      8----0----4---10     20-_  |  _16---22     32---24---28---34       
      |    |    |    |      |  -12_-  |    |      |    |    |    |       
      | 0  |  1 | 4  |      |    |    |    |      | 6  |  7 | 10 |
      1----2----3----9     13---14---15---21     25---26---27---33
   */
   inline void create_twelve_hex_patch(PatchData &pd, MsqError &err) 
   {
      pd.reserve_vertex_capacity(36, err); MSQ_CHKERR(err);
      pd.add_vertex(NULL, NULL, 1, 1, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 0, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 0, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 0, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 1, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 2, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 2, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 2, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 1, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 0, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 1, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 2, -1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      
      pd.add_vertex(NULL, NULL, 1,.5, 0, true, err, MsqVertex::MSQ_NO_VTX_FLAG); 
      pd.add_vertex(NULL, NULL, 0, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 1, 0, true, err, MsqVertex::MSQ_NO_VTX_FLAG); 
      pd.add_vertex(NULL, NULL, 2, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 1, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 0, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 1, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 2, 0, true, err, MsqVertex::MSQ_HARD_FIXED); 
      
      pd.add_vertex(NULL, NULL, 1, 1, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 0, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 0, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 0, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 1, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 2, 2, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 1, 2, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 2, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 0, 1, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 0, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 1, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      pd.add_vertex(NULL, NULL, 3, 2, 1, true, err, MsqVertex::MSQ_HARD_FIXED); 
      
      size_t ind[8];
      pd.reserve_element_capacity(12, err); MSQ_CHKERR(err);
      ind[0]=1; ind[1]=2; ind[2]=0; ind[3]=8; ind[4]=13; ind[5]=14; ind[6]=12; ind[7]=20; // 0
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=2; ind[1]=3; ind[2]=4; ind[3]=0; ind[4]=14; ind[5]=15; ind[6]=16; ind[7]=12; // 1
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=8; ind[1]=0; ind[2]=6; ind[3]=7; ind[4]=20; ind[5]=12; ind[6]=18; ind[7]=19; // 2
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=0; ind[1]=4; ind[2]=5; ind[3]=6; ind[4]=12; ind[5]=16; ind[6]=17; ind[7]=18; // 3
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=3; ind[1]=9; ind[2]=10; ind[3]=4; ind[4]=15; ind[5]=21; ind[6]=22; ind[7]=16; // 4
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=4; ind[1]=10; ind[2]=11; ind[3]=5; ind[4]=16; ind[5]=22; ind[6]=23; ind[7]=17; // 5
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=13; ind[1]=14; ind[2]=12; ind[3]=20; ind[4]=25; ind[5]=26; ind[6]=24; ind[7]=32; // 6
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=14; ind[1]=15; ind[2]=16; ind[3]=12; ind[4]=26; ind[5]=27; ind[6]=28; ind[7]=24; // 7
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=20; ind[1]=12; ind[2]=18; ind[3]=19; ind[4]=32; ind[5]=24; ind[6]=30; ind[7]=31; // 8
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=12; ind[1]=16; ind[2]=17; ind[3]=18; ind[4]=24; ind[5]=28; ind[6]=29; ind[7]=30; // 9
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=15; ind[1]=21; ind[2]=22; ind[3]=16; ind[4]=27; ind[5]=33; ind[6]=34; ind[7]=28; // 10
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
      ind[0]=16; ind[1]=22; ind[2]=23; ind[3]=17; ind[4]=28; ind[5]=34; ind[6]=35; ind[7]=29; // 11
      pd.add_element(NULL, NULL, ind, HEXAHEDRON, err); MSQ_CHKERR(err);
   }

   /* Patch used in several quality metric tests.
      Our triangular patch is made of two tris.  tri_1 is a perfect
      equilateral (the ideal for most metrics).  tri_2 is an arbitrary
      triangle.
   */
   inline void create_qm_two_tri_patch(PatchData &triPatch, MsqError &err) 
   {
     size_t ind[20];
     size_t elem_ind[8];
     triPatch.reserve_vertex_capacity(4, err); MSQ_CHKERR(err);
     ind[0]=triPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=triPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=triPatch.add_vertex(NULL, NULL, .5,sqrt(3.0)/2.0,0.0,
                                true, err);
     MSQ_CHKERR(err);
     ind[3]=triPatch.add_vertex(NULL, NULL, 2.0,2.0,2.0, true, err);
     MSQ_CHKERR(err);
     triPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal equilateral
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     triPatch.add_element(NULL, NULL, elem_ind, TRIANGLE, err);
     MSQ_CHKERR(err);
       //add "arbitrary" tri
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[3];
     elem_ind[2] = ind[1];
     triPatch.add_element(NULL, NULL, elem_ind, TRIANGLE, err);
     MSQ_CHKERR(err);

   }
     /* Patch used in several quality metric tests.
       Our quad patch is made of two quads.  quad_1 is a perfect
       square (the ideal for most metrics).  quad_2 is an arbitrary
       quad.
     */
   inline void create_qm_two_quad_patch(PatchData &quadPatch, MsqError &err)
   {
     size_t ind[20];
     size_t elem_ind[8];
     quadPatch.reserve_vertex_capacity(6, err); MSQ_CHKERR(err);
     ind[0]=quadPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=quadPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=quadPatch.add_vertex(NULL, NULL, 1.0,1.0,0.0,true, err);
     MSQ_CHKERR(err);
     ind[3]=quadPatch.add_vertex(NULL, NULL, 0.0,1.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[4]=quadPatch.add_vertex(NULL, NULL, 2.0,-1.0,.5, true, err);
     MSQ_CHKERR(err);
     ind[5]=quadPatch.add_vertex(NULL, NULL, 1.5,1.0,1.0, true, err);
     MSQ_CHKERR(err);
     quadPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal quad
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     
     quadPatch.add_element(NULL, NULL, elem_ind, QUADRILATERAL, err);
     MSQ_CHKERR(err);
       //add "arbitrary" quad
     elem_ind[0] = ind[1];
     elem_ind[1] = ind[4];
     elem_ind[2] = ind[5];
     elem_ind[3] = ind[2];
     
     quadPatch.add_element(NULL, NULL, elem_ind, QUADRILATERAL, err);
     MSQ_CHKERR(err);
   }
     /* Patch used in several quality metric tests.
        Our tet patch is made of two tets.  tet_1 is a perfect
        equilateral (the ideal for most metrics).  tet_2 is an arbitrary
        tet.
     */
   inline void create_qm_two_tet_patch(PatchData &tetPatch, MsqError &err)
   {
     size_t ind[20];
     size_t elem_ind[8];
     tetPatch.reserve_vertex_capacity(5, err); MSQ_CHKERR(err);
     ind[0]=tetPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=tetPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=tetPatch.add_vertex(NULL, NULL, 0.5,sqrt(3.0)/2.0,0.0,true, err);
     MSQ_CHKERR(err);
     ind[3]=tetPatch.add_vertex(NULL, NULL, .5, sqrt(3.0)/6.0,
                                 sqrt(2.0)/sqrt(3.0), true, err);
     MSQ_CHKERR(err);
     ind[4]=tetPatch.add_vertex(NULL, NULL, 2.0,3.0,-.5, true, err);
     MSQ_CHKERR(err);
    
     tetPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal tet
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     
     tetPatch.add_element(NULL, NULL, elem_ind, TETRAHEDRON, err);
     MSQ_CHKERR(err);
       //add "arbitrary" tet
     elem_ind[0] = ind[1];
     elem_ind[1] = ind[4];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     
     tetPatch.add_element(NULL, NULL, elem_ind, TETRAHEDRON, err);
     MSQ_CHKERR(err);
   }
   /* Patch used in seveal quality metric tests.
      Our hex patch is made of two hexes.  hex_1 is a perfect
      unit cube (the ideal for most metrics).  hex_2 is an arbitrary
      hex.
   */
   inline void create_qm_two_hex_patch(PatchData &hexPatch, MsqError &err)
   {
     size_t ind[20];
     size_t elem_ind[8];
     hexPatch.reserve_vertex_capacity(12, err); MSQ_CHKERR(err);
     ind[0]=hexPatch.add_vertex(NULL, NULL, 0.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[1]=hexPatch.add_vertex(NULL, NULL, 1.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[2]=hexPatch.add_vertex(NULL, NULL, 1.0, 1.0, 0.0, true, err);
     MSQ_CHKERR(err);
     ind[3]=hexPatch.add_vertex(NULL, NULL, 0.0, 1.0, 0.0, true, err);
     MSQ_CHKERR(err);
     ind[4]=hexPatch.add_vertex(NULL, NULL, 0.0,0.0,1.0, true, err);
     MSQ_CHKERR(err);
     ind[5]=hexPatch.add_vertex(NULL, NULL, 1.0,0.0,1.0, true, err);
     MSQ_CHKERR(err);
     ind[6]=hexPatch.add_vertex(NULL, NULL, 1.0, 1.0, 1.0, true, err);
     MSQ_CHKERR(err);
     ind[7]=hexPatch.add_vertex(NULL, NULL, 0.0, 1.0, 1.0, true, err);
     MSQ_CHKERR(err);
     ind[8]=hexPatch.add_vertex(NULL, NULL, 2.0,0.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[9]=hexPatch.add_vertex(NULL, NULL, 2.0,1.0,0.0, true, err);
     MSQ_CHKERR(err);
     ind[10]=hexPatch.add_vertex(NULL, NULL, 2.0, -1.0, 1.0, true, err);
     MSQ_CHKERR(err);
     ind[11]=hexPatch.add_vertex(NULL, NULL, 3.0, 2.0, 1.0, true, err);
     MSQ_CHKERR(err);
    
     hexPatch.reserve_element_capacity(2, err); MSQ_CHKERR(err);
     MSQ_CHKERR(err);
       //add ideal hex
     elem_ind[0] = ind[0];
     elem_ind[1] = ind[1];
     elem_ind[2] = ind[2];
     elem_ind[3] = ind[3];
     elem_ind[4] = ind[4];
     elem_ind[5] = ind[5];
     elem_ind[6] = ind[6];
     elem_ind[7] = ind[7];
     
     hexPatch.add_element(NULL, NULL, elem_ind, HEXAHEDRON, err);
     MSQ_CHKERR(err);
       //add "arbitrary" hex
     elem_ind[0] = ind[1];
     elem_ind[1] = ind[8];
     elem_ind[2] = ind[9];
     elem_ind[3] = ind[2];
     elem_ind[4] = ind[5];
     elem_ind[5] = ind[10];
     elem_ind[6] = ind[11];
     elem_ind[7] = ind[6];
     
     hexPatch.add_element(NULL, NULL, elem_ind, HEXAHEDRON, err);
     MSQ_CHKERR(err);
   }
   
} // namespace

#endif // PatchDataInstances_hpp
