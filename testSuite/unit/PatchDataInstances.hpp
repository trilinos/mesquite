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
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 12-Nov-02 at 18:05:56
//  LAST-MOD:  9-Jun-04 at 14:43:39 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PatchDataInstances.hpp

This header file contains some functions to instantiates particular PatchData Objects.
Those objects can be used in unit tests.
Patches must be allocated and dealocated by the caller. 

\author Thomas Leurent
\author Michael Brewer

*/
// DESCRIP-END.
//

#ifndef PatchDataInstances_hpp
#define PatchDataInstances_hpp

#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MeshSet.hpp"
#include "PlanarDomain.hpp"

#include <math.h>
#include <iostream>

#include "cppunit/extensions/HelperMacros.h"

namespace Mesquite
{
  //! must be called in sync with create_...._patch_with_domain
  inline void destroy_patch_with_domain(PatchData &pd)
  {
    delete pd.get_mesh_set()->get_domain_constraint();
    delete pd.get_mesh_set();
  }


  /*! creates a patch containing one ideal hexahedra
  */
   inline void create_one_hex_patch(PatchData &one_hex_patch, MsqError &err)
   {
       // creates empty Patch
     one_hex_patch.set_num_vertices(8);
     
       // Fills up with vertices for ideal hexahedra.
     double coords[3];
     MsqVertex* vert_array = one_hex_patch.get_vertex_array(err);
     
     coords[0] = 1.0; coords[1] = 1.0; coords[2] = 1.0;
     vert_array[0] = coords;
     
     coords[0] = 2; coords[1] = 1; coords[2] = 1;
     vert_array[1] = coords;
     
     coords[0] = 2.; coords[1] = 2.; coords[2] = 1;
     vert_array[2] = coords;
     
     coords[0] = 1.; coords[1] = 2.; coords[2] = 1;
     vert_array[3] = coords;
     
     coords[0] = 1.; coords[1] = 1.; coords[2] = 2;
     vert_array[4] = coords;
     
     coords[0] = 2.; coords[1] = 1.; coords[2] = 2;
     vert_array[5] = coords;
     
     coords[0] = 2.; coords[1] = 2.; coords[2] = 2;
     vert_array[6] = coords;
     
     coords[0] = 1.; coords[1] = 2.; coords[2] = 2;
     vert_array[7] = coords;
     
       // patch has only one element: an ideal hex.
     one_hex_patch.set_num_elements(1);
     size_t indices[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
     MsqMeshEntity& hex = one_hex_patch.element_by_index(0);
     hex.set_element_type(Mesquite::HEXAHEDRON);
     memcpy(hex.get_modifiable_vertex_index_array(),
            indices,
            8*sizeof(size_t));
   }
   

   //! creates a Patch containing an ideal tetrahedra
   inline void create_one_tet_patch(PatchData &one_tet_patch, MsqError &err)
   {
       /* *********************FILL TET************************ */
       // creates empty Patch
     one_tet_patch.set_num_vertices(4);
     one_tet_patch.set_num_elements(1);
     
       // Fills up with vertices for ideal tet
     double coords[3];
     MsqVertex* vert_array = one_tet_patch.get_vertex_array(err);
     
     coords[0] = 1; coords[1] = 1; coords[2] = 1;
     vert_array[0] = coords;
     
     coords[0] = 2; coords[1] = 1; coords[2] = 1;
     vert_array[1] = coords;

     coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/2.0; coords[2] = 1;
     vert_array[2] = coords;

     coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/6.0;
     coords[2] = 1+sqrt(2.0)/sqrt(3.0);
     vert_array[3] = coords;
     
       // patch has only one element: an ideal tet.
     size_t indices[4] = { 0, 1, 2, 3 };
     MsqMeshEntity& tet = one_tet_patch.element_by_index(0);
     tet.set_element_type(Mesquite::TETRAHEDRON);
     memcpy(tet.get_modifiable_vertex_index_array(),
            indices,
            4*sizeof(size_t));
   }
      //! creates a Patch containing an ideal tetrahedra, inverted
   inline void create_one_inverted_tet_patch(PatchData &one_tet_patch,
                                             MsqError &err)
   {
       /* *********************FILL TET************************* */
       // creates empty Patch
     one_tet_patch.set_num_vertices(4);
     one_tet_patch.set_num_elements(1);
     
       // Fills up with vertices for ideal tet
     double coords[3];
     MsqVertex* vert_array = one_tet_patch.get_vertex_array(err);
     
     coords[0] = 1; coords[1] = 1; coords[2] = 1;
     vert_array[0] = coords;
     
     coords[0] = 2; coords[1] = 1; coords[2] = 1;
     vert_array[1] = coords;

     coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/2.0; coords[2] = 1;
     vert_array[2] = coords;

     coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/6.0;
     coords[2] = 1-sqrt(2.0)/sqrt(3.0);
     vert_array[3] = coords;
     
       // patch has only one element: an ideal tet.
     size_t indices[4] = { 0, 1, 2, 3 };
     MsqMeshEntity& tet = one_tet_patch.element_by_index(0);
     tet.set_element_type(Mesquite::TETRAHEDRON);
     memcpy(tet.get_modifiable_vertex_index_array(),
            indices,
            4*sizeof(size_t));
   }
   
   //! creates a Patch containing an ideal quadrilateral
   inline void create_one_quad_patch(PatchData &one_qua_patch, MsqError &err)
   {
       /* *********************FILL QUAD************************* */
     one_qua_patch.set_num_vertices(4);
     one_qua_patch.set_num_elements(1);
     
       // Fills up with vertices for ideal quad
     double coords[3];
     MsqVertex* vert_array = one_qua_patch.get_vertex_array(err);
     
     coords[0] = 1; coords[1] = 1; coords[2] = 1;
     vert_array[0] = coords;
     coords[0] = 2; coords[1] = 1; coords[2] = 1;
     vert_array[1] = coords;
     coords[0] = 2; coords[1] = 2; coords[2] = 1;
     vert_array[2] = coords;
     coords[0] = 1; coords[1] = 2 ; coords[2] = 1;
     vert_array[3] = coords;
     
       // patch has only one element: an ideal quad.
     size_t indices[4] = { 0, 1, 2, 3 };
     MsqMeshEntity& quad = one_qua_patch.element_by_index(0);
     quad.set_element_type(Mesquite::QUADRILATERAL);
     memcpy(quad.get_modifiable_vertex_index_array(),
            indices,
            4*sizeof(size_t));
   }
   
   
     /*! \fn create_one_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing an ideal triangle
          /   \
         0-----1
         This Patch also has the normal information. 
     */
   inline void create_one_tri_patch(PatchData &one_tri_patch, MsqError &err)
   {
       /* ************** Creates normal info ******************* */
     MeshSet* mesh_set1 = new MeshSet;;
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt, NULL);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     one_tri_patch.set_mesh_set(mesh_set1);

       /* *********************FILL tri************************* */
     one_tri_patch.set_num_vertices(3);
     one_tri_patch.set_num_elements(1);
     
       // Fills up with vertices for ideal tri
     double coords[3];
     MsqVertex* vert_array = one_tri_patch.get_vertex_array(err);
     
     coords[0] = 1; coords[1] = 1; coords[2] = 1;
     vert_array[0] = coords;
     
     coords[0] = 2; coords[1] = 1; coords[2] = 1;
     vert_array[1] = coords;
     
     coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/2.0; coords[2] = 1;
     vert_array[2] = coords;
     
       // patch has only one element: an ideal tri
     size_t indices[3] = { 0, 1, 2 };
     MsqMeshEntity& tri = one_tri_patch.element_by_index(0);
     tri.set_element_type(Mesquite::TRIANGLE);
     memcpy(tri.get_modifiable_vertex_index_array(),
            indices,
            3*sizeof(size_t));
   }
     
//     inline void create_one_inverted_tri_patch(PatchData &one_tri_patch, MsqError &err)
//    {

//      // **********************FILL tri*************************
//        // creates empty Patch
//      one_tri_patch.set_num_vertices(3);
//      one_tri_patch.set_num_elements(1);
     
//        // Fills up with vertices for ideal tri
//      double coords[3];
//      MsqVertex* vert_array = one_tri_patch.get_vertex_array(err);
     
//      coords[0] = 1; coords[1] = 1; coords[2] = 1;
//      vert_array[0] = coords;
     
//      coords[0] = 2; coords[1] = 1; coords[2] = 1;
//      vert_array[1] = coords;
     
//      coords[0] = 1.5; coords[1] = 1-sqrt(3.0)/2.0; coords[2] = 1;
//      vert_array[2] = coords;
     
//        // patch has only one element: an ideal tri
//      size_t indices[3] = { 0, 1, 2 };
//      MsqMeshEntity& tri = one_tri_patch.element_by_index(0);
//      tri.set_element_type(Mesquite::TRIANGLE);
//      memcpy(tri.get_modifiable_vertex_index_array(),
//             indices,
//             3*sizeof(size_t));
//    }
      
  
   /*! \fn create_two_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing two ideal triangles
          / 0 \
         0-----1
          \ 1 /
           \ /
            3
         This Patch also has the normal information. 
   */
   inline void create_two_tri_patch(PatchData &pd, MsqError &err)
   {
       /* ************** Creates normal info ******************* */
     MeshSet* mesh_set1 = new MeshSet;;
     Vector3D pnt(0,0,1);
     Vector3D s_norm(0,0,3);
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt, NULL);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     pd.set_mesh_set(mesh_set1);

       // **********************FILL tri*************************
       // creates empty Patch
     pd.set_num_vertices(4);
     pd.set_num_elements(2);
     
       // Fills up with vertices for ideal triangles
     double coords[3];
     MsqVertex* vert_array = pd.get_vertex_array(err);
     
     coords[0] = 1; coords[1] = 1; coords[2] = 1;
     vert_array[0].set(coords);
     
     coords[0] = 2; coords[1] = 1; coords[2] = 1;
     vert_array[1].set(coords);
     
     coords[0] = 1.5; coords[1] = 1+sqrt(3.0)/2.0; coords[2] = 1;
     vert_array[2].set(coords);
     
     coords[0] = 1.5; coords[1] = 1-sqrt(3.0)/2.0; coords[2] = 1;
     vert_array[3].set(coords);
    
     size_t indices[3] = { 0, 1, 2 };
     MsqMeshEntity& tri = pd.element_by_index(0);
     tri.set_element_type(Mesquite::TRIANGLE);
     memcpy(tri.get_modifiable_vertex_index_array(),
            indices,
            3*sizeof(size_t));
     
     indices[0] = 0; indices[1] = 3; indices[2] = 1;
     MsqMeshEntity& tri2 = pd.element_by_index(1);
     tri2.set_element_type(Mesquite::TRIANGLE);
     memcpy(tri2.get_modifiable_vertex_index_array(),
            indices,
            3*sizeof(size_t));
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
     four_quads.set_num_vertices(9);
     four_quads.set_num_elements(4);
     
       // Fill up vertices
     double coords[3];
     MsqVertex* vert_array = four_quads.get_vertex_array(err);
     
     coords[0] = 1; coords[1] = .5, coords[2] = 0;
     vert_array[0] = coords;
     
     coords[0] = 0, coords[1] = 0, coords[2] = 0;
     vert_array[1] = coords;
     
     coords[0] = 1, coords[1] = 0, coords[2] = 0;
     vert_array[2] = coords;
     
     coords[0] = 2, coords[1] = 0, coords[2] = 0;
     vert_array[3] = coords;
     
     coords[0] = 2, coords[1] = 1, coords[2] = 0;
     vert_array[4] = coords;
     
     coords[0] = 2, coords[1] = 2, coords[2] = 0;
     vert_array[5] = coords;
     
     coords[0] = 1, coords[1] = 2, coords[2] = 0;
     vert_array[6] = coords;
     
     coords[0] = 0, coords[1] = 2, coords[2] = 0;
     vert_array[7] = coords;
     
     coords[0] = 0, coords[1] = 1, coords[2] = 0;
     vert_array[8] = coords;
      
     size_t ind[4];
     ind[0] = 1; ind[1]=2; ind[2]=0; ind[3]=8;
     four_quads.element_by_index(0).set_element_type(Mesquite::QUADRILATERAL);
     memcpy(four_quads.element_by_index(0).get_modifiable_vertex_index_array(),
            ind,
            4*sizeof(size_t));
     
     ind[0] = 2; ind[1]=3; ind[2]=4; ind[3]=0;
     four_quads.element_by_index(1).set_element_type(Mesquite::QUADRILATERAL);
     memcpy(four_quads.element_by_index(1).get_modifiable_vertex_index_array(),
            ind,
            4*sizeof(size_t));
     
     ind[0] = 8; ind[1]=0; ind[2]=6; ind[3]=7;
     four_quads.element_by_index(2).set_element_type(Mesquite::QUADRILATERAL);
     memcpy(four_quads.element_by_index(2).get_modifiable_vertex_index_array(),
            ind,
            4*sizeof(size_t));
     
     ind[0] = 0; ind[1]=4; ind[2]=5; ind[3]=6;
     four_quads.element_by_index(3).set_element_type(Mesquite::QUADRILATERAL);
     memcpy(four_quads.element_by_index(3).get_modifiable_vertex_index_array(),
            ind,
            4*sizeof(size_t));
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

      use destroy_patch_with_domain() in sync.
   */
   inline void create_six_quads_patch_with_domain(PatchData &pd, MsqError &err) 
   {
     // associates domain
     MeshSet* mesh_set1 = new MeshSet;;
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt, NULL);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     pd.set_mesh_set(mesh_set1);

     pd.set_num_vertices(12);
     MsqVertex* vert_array = pd.get_vertex_array(err);
     vert_array[0] = Vector3D(1,.5, 0);
     vert_array[1] = Vector3D(0, 0, 0);
     vert_array[1].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[2] = Vector3D(1, 0, 0);
     vert_array[2].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[3] = Vector3D(2, 0, 0);
     vert_array[3].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[4] = Vector3D(2, 1, 0);
     vert_array[5] = Vector3D(2, 2, 0);
     vert_array[5].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[6] = Vector3D(1, 2, 0);
     vert_array[6].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[7] = Vector3D(0, 2, 0);
     vert_array[7].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[8] = Vector3D(0, 1, 0);
     vert_array[8].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[9] = Vector3D(3, 0, 0);
     vert_array[9].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[10] = Vector3D(3, 1, 0);
     vert_array[10].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[11] = Vector3D(3, 2, 0);
     vert_array[11].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     
     size_t ind[6][4] =
     { { 1,  2,  0, 8},
       { 2,  3,  4, 0},
       { 8,  0,  6, 7},
       { 0,  4,  5, 6},
       { 3,  9, 10, 4},
       { 4, 10, 11, 5}
     };
     
     pd.set_num_elements(6);
     for (int i = 0; i < 6; i++)
     {
       pd.element_by_index(i).set_element_type(QUADRILATERAL);
       size_t *ind_array =
          pd.element_by_index(i).get_modifiable_vertex_index_array();
       for (int j = 0; j < 4; j++)
         ind_array[j] = ind[i][j];
     }
   }
   

   /*! \fn create_six_quads_patch_inverted_with_domain(PatchData &four_quads, MsqError &err)
     our 2D set up: 6 quads, 1 center vertex outcentered by (0,-0.5), the other centered
      7____6____5___11
      |    |    |    |
      | 2  |  3 | 5  |
      8    |    4---10       vertex 1 is at (0,0)
      |\       /|    |       vertex 11 is at (3,2)
      |    |    | 4  |
      1----2----3----9
         \  /
          0      
      use destroy_patch_with_domain() in sync.
   */
   inline void create_six_quads_patch_inverted_with_domain(PatchData &pd, MsqError &err) 
   {
     create_six_quads_patch_with_domain(pd,err); MSQ_CHKERR(err);

     MsqVertex* vtx = pd.get_vertex_array(err); 

     Vector3D displacement(0,-1.5,0);
     
     vtx[0] += displacement;
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
     pd.set_num_vertices(36);
     MsqVertex* verts = pd.get_vertex_array(err); MSQ_CHKERR(err);
     
     verts[0].set(1, 1, -1);
     verts[0].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[1].set(0, 0, -1);
     verts[1].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[2].set(1, 0, -1);
     verts[2].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[3].set(2, 0, -1);
     verts[3].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[4].set(2, 1, -1);
     verts[4].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[5].set(2, 2, -1);
     verts[5].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[6].set(1, 2, -1);
     verts[6].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[7].set(0, 2, -1);
     verts[7].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[8].set(0, 1, -1);
     verts[8].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[9].set(3, 0, -1);
     verts[9].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[10].set(3, 1, -1);
     verts[10].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[11].set(3, 2, -1);
     verts[11].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
      
     verts[12].set(1,.5, 0);
     verts[12].set_vertex_flag(MsqVertex::MSQ_NO_VTX_FLAG); 
     verts[13].set(0, 0, 0);
     verts[13].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[14].set(1, 0, 0);
     verts[14].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[15].set(2, 0, 0);
     verts[15].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[16].set(2, 1, 0);
     verts[16].set_vertex_flag(MsqVertex::MSQ_NO_VTX_FLAG); 
     verts[17].set(2, 2, 0);
     verts[17].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[18].set(1, 2, 0);
     verts[18].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[19].set(0, 2, 0);
     verts[19].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[20].set(0, 1, 0);
     verts[20].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[21].set(3, 0, 0);
     verts[21].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[22].set(3, 1, 0);
     verts[22].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[23].set(3, 2, 0);
     verts[23].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
      
     verts[24].set(1, 1, 1);
     verts[24].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[25].set(0, 0, 1);
     verts[25].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[26].set(1, 0, 1);
     verts[26].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[27].set(2, 0, 1);
     verts[27].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[28].set(2, 1, 1);
     verts[28].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[29].set(2, 2, 1);
     verts[29].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[30].set(1, 2, 1);
     verts[30].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[31].set(0, 2, 1);
     verts[31].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[32].set(0, 1, 1);
     verts[32].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[33].set(3, 0, 1);
     verts[33].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[34].set(3, 1, 1);
     verts[34].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[35].set(3, 2, 1);
     verts[35].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     
     size_t ind[8];
     pd.set_num_elements(12);
     ind[0]=1; ind[1]=2; ind[2]=0; ind[3]=8; ind[4]=13; ind[5]=14; ind[6]=12; ind[7]=20; // 0
     pd.element_by_index(0).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(0).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=2; ind[1]=3; ind[2]=4; ind[3]=0; ind[4]=14; ind[5]=15; ind[6]=16; ind[7]=12; // 1
     pd.element_by_index(1).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(1).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=8; ind[1]=0; ind[2]=6; ind[3]=7; ind[4]=20; ind[5]=12; ind[6]=18; ind[7]=19; // 2
     pd.element_by_index(2).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(2).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=0; ind[1]=4; ind[2]=5; ind[3]=6; ind[4]=12; ind[5]=16; ind[6]=17; ind[7]=18; // 3
     pd.element_by_index(3).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(3).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=3; ind[1]=9; ind[2]=10; ind[3]=4; ind[4]=15; ind[5]=21; ind[6]=22; ind[7]=16; // 4
     pd.element_by_index(4).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(4).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=4; ind[1]=10; ind[2]=11; ind[3]=5; ind[4]=16; ind[5]=22; ind[6]=23; ind[7]=17; // 5
     pd.element_by_index(5).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(5).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=13; ind[1]=14; ind[2]=12; ind[3]=20; ind[4]=25; ind[5]=26; ind[6]=24; ind[7]=32; // 6
     pd.element_by_index(6).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(6).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=14; ind[1]=15; ind[2]=16; ind[3]=12; ind[4]=26; ind[5]=27; ind[6]=28; ind[7]=24; // 7
     pd.element_by_index(7).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(7).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=20; ind[1]=12; ind[2]=18; ind[3]=19; ind[4]=32; ind[5]=24; ind[6]=30; ind[7]=31; // 8
     pd.element_by_index(8).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(8).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=12; ind[1]=16; ind[2]=17; ind[3]=18; ind[4]=24; ind[5]=28; ind[6]=29; ind[7]=30; // 9
     pd.element_by_index(9).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(9).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=15; ind[1]=21; ind[2]=22; ind[3]=16; ind[4]=27; ind[5]=33; ind[6]=34; ind[7]=28; // 10
     pd.element_by_index(10).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(10).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
     ind[0]=16; ind[1]=22; ind[2]=23; ind[3]=17; ind[4]=28; ind[5]=34; ind[6]=35; ind[7]=29; // 11
     pd.element_by_index(11).set_element_type(HEXAHEDRON);
     memcpy(pd.element_by_index(11).get_modifiable_vertex_index_array(),
            ind,
            8*sizeof(size_t));
   }
   
   inline void create_twelve_hex_patch_inverted(PatchData &pd, MsqError &err)
   {
     create_twelve_hex_patch(pd,err); MSQ_CHKERR(err); 

     MsqVertex* vtx = pd.get_vertex_array(err); 

     Vector3D displacement(0,0,1.5);
     
     vtx[16] += displacement;
   }
     
     
   /* Patch used in several quality metric tests.
      Our triangular patch is made of two tris.  tri_1 is a perfect
      equilateral (the ideal for most metrics).  tri_2 is an arbitrary
      triangle.
      Memory allocated in this function must be deallocated with
      destroy_patch_with_domain().
   */
   inline void create_qm_two_tri_patch_with_domain(PatchData &triPatch, MsqError &err) 
   {
     MeshSet* mesh_set1 = new MeshSet;;
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt, NULL);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     triPatch.set_mesh_set(mesh_set1);

     triPatch.set_num_vertices(4);
     triPatch.vertex_by_index(0).set(0.0, 0.0, 0.0);
     triPatch.vertex_by_index(1).set(1.0, 0.0, 0.0);
     triPatch.vertex_by_index(2).set(.5, sqrt(3.0)/2.0, 0.0);
     triPatch.vertex_by_index(3).set(2.0, -4.0, 2.0);
     
     triPatch.set_num_elements(2);
     
       //add ideal equilateral
     triPatch.element_by_index(0).set_element_type(TRIANGLE);
     triPatch.element_by_index(0).set_vertex_index(0, 0);
     triPatch.element_by_index(0).set_vertex_index(1, 1);
     triPatch.element_by_index(0).set_vertex_index(2, 2);
     
       //add "arbitrary" tri
     triPatch.element_by_index(1).set_element_type(TRIANGLE);
     triPatch.element_by_index(1).set_vertex_index(0, 0);
     triPatch.element_by_index(1).set_vertex_index(1, 3);
     triPatch.element_by_index(1).set_vertex_index(2, 1);
   }
   
     /* Patch used in several quality metric tests.
       Our quad patch is made of two quads.  quad_1 is a perfect
       square (the ideal for most metrics).  quad_2 is an arbitrary
       quad.
       Memory allocated in this function must be deallocated with
       destroy_patch_with_domain().
     */
   inline void create_qm_two_quad_patch_with_domain(PatchData &quadPatch, MsqError &err)
   {
     MeshSet* mesh_set1 = new MeshSet;;
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt, NULL);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     quadPatch.set_mesh_set(mesh_set1);

     // Add vertices
     quadPatch.set_num_vertices(6);
     quadPatch.vertex_by_index(0).set(0.0, 0.0, 0.0);
     quadPatch.vertex_by_index(1).set(1.0, 0.0, 0.0);
     quadPatch.vertex_by_index(2).set(1.0, 1.0, 0.0);
     quadPatch.vertex_by_index(3).set(0.0, 1.0, 0.0);
     quadPatch.vertex_by_index(4).set(2.0, -1.0, .5);
     quadPatch.vertex_by_index(5).set(1.5, 1.0, 1.0);
     
     quadPatch.set_num_elements(2);
     
       //add ideal quad
     quadPatch.element_by_index(0).set_element_type(QUADRILATERAL);
     quadPatch.element_by_index(0).set_vertex_index(0, 0);
     quadPatch.element_by_index(0).set_vertex_index(1, 1);
     quadPatch.element_by_index(0).set_vertex_index(2, 2);
     quadPatch.element_by_index(0).set_vertex_index(3, 3);
     
       //add "arbitrary" quad
     quadPatch.element_by_index(1).set_element_type(QUADRILATERAL);
     quadPatch.element_by_index(1).set_vertex_index(0, 1);
     quadPatch.element_by_index(1).set_vertex_index(1, 4);
     quadPatch.element_by_index(1).set_vertex_index(2, 5);
     quadPatch.element_by_index(1).set_vertex_index(3, 2);
   }
  
     /* Patch used in several quality metric tests.
        Our tet patch is made of two tets.  tet_1 is a perfect
        equilateral (the ideal for most metrics).  tet_2 is an arbitrary
        tet.
     */
   inline void create_qm_two_tet_patch(PatchData &tetPatch, MsqError &/*err*/)
   {
     tetPatch.set_num_vertices(5);

     tetPatch.vertex_by_index(0).set(0.0, 0.0, 0.0);
     tetPatch.vertex_by_index(1).set(1.0, 0.0, 0.0);
     tetPatch.vertex_by_index(2).set(0.5, sqrt(3.0)/2.0, 0.0);
     tetPatch.vertex_by_index(3).set(.5, sqrt(3.0)/6.0, sqrt(2.0)/sqrt(3.0));
     tetPatch.vertex_by_index(4).set(2.0, 3.0, -.5);
     
     tetPatch.set_num_elements(2);
     
       //add ideal tet
     tetPatch.element_by_index(0).set_element_type(TETRAHEDRON);
     tetPatch.element_by_index(0).set_vertex_index(0, 0);
     tetPatch.element_by_index(0).set_vertex_index(1, 1);
     tetPatch.element_by_index(0).set_vertex_index(2, 2);
     tetPatch.element_by_index(0).set_vertex_index(3, 3);
     
       //add "arbitrary" tet
     tetPatch.element_by_index(1).set_element_type(TETRAHEDRON);
     tetPatch.element_by_index(1).set_vertex_index(0, 1);
     tetPatch.element_by_index(1).set_vertex_index(1, 4);
     tetPatch.element_by_index(1).set_vertex_index(2, 2);
     tetPatch.element_by_index(1).set_vertex_index(3, 3);
   }
   /* Patch used in seveal quality metric tests.
      Our hex patch is made of two hexes.  hex_1 is a perfect
      unit cube (the ideal for most metrics).  hex_2 is an arbitrary
      hex.
   */
   inline void create_qm_two_hex_patch(PatchData &hexPatch, MsqError &/*err*/)
   {
     hexPatch.set_num_vertices(12);

     hexPatch.vertex_by_index(0).set(0.0,0.0,0.0);
     hexPatch.vertex_by_index(1).set(1.0,0.0,0.0);
     hexPatch.vertex_by_index(2).set(1.0, 1.0, 0.0);
     hexPatch.vertex_by_index(3).set(0.0, 1.0, 0.0);
     hexPatch.vertex_by_index(4).set(0.0,0.0,1.0);
     hexPatch.vertex_by_index(5).set(1.0,0.0,1.0);
     hexPatch.vertex_by_index(6).set(1.0, 1.0, 1.0);
     hexPatch.vertex_by_index(7).set(0.0, 1.0, 1.0);
     hexPatch.vertex_by_index(8).set(2.0,0.0,0.0);
     hexPatch.vertex_by_index(9).set(2.0,1.0,0.0);
     hexPatch.vertex_by_index(10).set(2.0, -1.0, 1.0);
     hexPatch.vertex_by_index(11).set(3.0, 2.0, 1.0);
     
     hexPatch.set_num_elements(2);
     
       //add ideal hex
     hexPatch.element_by_index(0).set_element_type(HEXAHEDRON);
     hexPatch.element_by_index(0).set_vertex_index(0, 0);
     hexPatch.element_by_index(0).set_vertex_index(1, 1);
     hexPatch.element_by_index(0).set_vertex_index(2, 2);
     hexPatch.element_by_index(0).set_vertex_index(3, 3);
     hexPatch.element_by_index(0).set_vertex_index(4, 4);
     hexPatch.element_by_index(0).set_vertex_index(5, 5);
     hexPatch.element_by_index(0).set_vertex_index(6, 6);
     hexPatch.element_by_index(0).set_vertex_index(7, 7);
     
       //add "arbitrary" hex
     hexPatch.element_by_index(1).set_element_type(HEXAHEDRON);
     hexPatch.element_by_index(1).set_vertex_index(0, 1);
     hexPatch.element_by_index(1).set_vertex_index(1, 8);
     hexPatch.element_by_index(1).set_vertex_index(2, 9);
     hexPatch.element_by_index(1).set_vertex_index(3, 2);
     hexPatch.element_by_index(1).set_vertex_index(4, 5);
     hexPatch.element_by_index(1).set_vertex_index(5, 10);
     hexPatch.element_by_index(1).set_vertex_index(6, 11);
     hexPatch.element_by_index(1).set_vertex_index(7, 6);
   }
   
} // namespace

#endif // PatchDataInstances_hpp
