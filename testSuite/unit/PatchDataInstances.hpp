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
//  LAST-MOD:  5-Dec-02 at 17:15:16 by Thomas Leurent
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
      int indices[8];
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
      int indices_tet[4];
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
      int indices_qua[4];
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
      int indices_tri[3];
      indices_tri[0] = 0; indices_tri[1] = 1; indices_tri[2] = 2;
      one_tri_patch.add_element(NULL, NULL, indices_tri, TRIANGLE, err);
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
      
      int ind[4];
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
     our 2D set up: 4 quads, center vertex outcentered by (0,-0.5)
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
      pd.reserve_vertex_capacity(9, err); MSQ_CHKERR(err);
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
      
      int ind[4];
      pd.reserve_element_capacity(6, err); MSQ_CHKERR(err);
      ind[0] = 1; ind[1]=2; ind[2]=0; ind[3]=8;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 2; ind[1]=3; ind[2]=4; ind[3]=0;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 8; ind[1]=0; ind[2]=6; ind[3]=7;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 0; ind[1]=4; ind[2]=5; ind[3]=6;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 3; ind[1]=9; ind[2]=10 ind[3]=4;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
      ind[0] = 4; ind[1]=10 ind[2]=11 ind[3]=5;
      pd.add_element(NULL, NULL, ind, QUADRILATERAL, err); MSQ_CHKERR(err);
   }


} // namespace

#endif // PatchDataInstances_hpp
