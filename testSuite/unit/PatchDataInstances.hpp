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
    if(pd.get_mesh_set()){
      if( pd.get_mesh_set()->get_domain_constraint() )
        delete pd.get_mesh_set()->get_domain_constraint();
      delete pd.get_mesh_set();
    }
    
  }
  
  inline void create_patch_mesh( PatchData& pd,
                                 size_t num_verts, 
                                 const double* coordinates,
                                 size_t num_elems,
                                 const size_t* connectivity,
                                 EntityTopology elem_type,
                                 MsqError& err )
  {
    size_t i;
    
    const size_t verts_per_elem = TopologyInfo::corners( elem_type );
    const size_t num_vert_uses = verts_per_elem * num_elems;
    pd.allocate_storage( num_verts, num_elems, num_vert_uses, err );
    MSQ_ERRRTN(err);
    
    for (i = 0; i < num_verts; ++i)
      pd.vertex_by_index(i) = coordinates + 3*i;
      
    for (i = 0; i < num_elems; ++i)
      pd.element_by_index(i).set_element_type( elem_type );
      
    memcpy( pd.get_connectivity_array(), connectivity, num_vert_uses * sizeof(size_t) );
    
    size_t* offsets = new size_t[num_elems+1];
    for (i = 0; i <= num_elems; ++i)
      offsets[i] = i * verts_per_elem;
    
    pd.initialize_data( offsets, err ); MSQ_CHKERR(err);
    delete [] offsets;
  }
    
    


  /*! creates a patch containing one ideal hexahedra
  */
   inline void create_one_hex_patch(PatchData &one_hex_patch, MsqError &err)
   {
     double coords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         2.0, 2.0, 1.0,
                         1.0, 2.0, 1.0,
                         1.0, 1.0, 2.0,
                         2.0, 1.0, 2.0,
                         2.0, 2.0, 2.0,
                         1.0, 2.0, 2.0 };
     
     size_t indices[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
     
     create_patch_mesh( one_hex_patch, 8, coords, 1, indices, HEXAHEDRON, err );
   }
   

   //! creates a Patch containing an ideal tetrahedra
   inline void create_one_tet_patch(PatchData &one_tet_patch, MsqError &err)
   {
     double coords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         1.5, 1+sqrt(3.0)/2.0, 1.0,
                         1.5, 1+sqrt(3.0)/6.0, 1+sqrt(2.0)/sqrt(3.0) };

     size_t indices[4] = { 0, 1, 2, 3 };

     create_patch_mesh( one_tet_patch, 4, coords, 1, indices, TETRAHEDRON, err );
   }

      //! creates a Patch containing an ideal tetrahedra, inverted
   inline void create_one_inverted_tet_patch(PatchData &one_tet_patch,
                                             MsqError &err)
   {
     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         1.5, 1+sqrt(3.0)/2.0, 1,
                         1.5, 1+sqrt(3.0)/6.0, 1-sqrt(2.0)/sqrt(3.0), };
                         
     size_t indices[4] = { 0, 1, 2, 3 };

     create_patch_mesh( one_tet_patch, 4, coords, 1, indices, TETRAHEDRON, err );
   }
   
   //! creates a Patch containing an ideal quadrilateral
   inline void create_one_quad_patch(PatchData &one_qua_patch, MsqError &err)
   {
     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         2, 2, 1,
                         1, 2 , 1 };
     
     size_t indices[4] = { 0, 1, 2, 3 };
     
     create_patch_mesh( one_qua_patch, 4, coords, 1, indices, QUADRILATERAL, err );
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
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     one_tri_patch.set_mesh_set(mesh_set1);

       /* *********************FILL tri************************* */
     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         1.5, 1+sqrt(3.0)/2.0, 1 };
     
     size_t indices[3] = { 0, 1, 2 };
     create_patch_mesh( one_tri_patch, 3, coords, 1, indices, TRIANGLE, err );
   }
     
  
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
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     pd.set_mesh_set(mesh_set1);

       // **********************FILL tri*************************

     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         1.5, 1+sqrt(3.0)/2.0, 1,
                         1.5, 1-sqrt(3.0)/2.0, 1 };

     size_t indices[] = { 0, 1, 2, 0, 3, 1 };
     
     create_patch_mesh( pd, 4, coords, 2, indices, TRIANGLE, err );
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
     double coords[] = { 1, .5, 0,
                         0, 0, 0,
                         1, 0, 0,
                         2, 0, 0,
                         2, 1, 0,
                         2, 2, 0,
                         1, 2, 0,
                         0, 2, 0,
                         0, 1, 0 };

     size_t indices[] = { 1, 2, 0, 8, 
                          2, 3, 4, 0,
                          8, 0, 6, 7,
                          0, 4, 5, 6 };
                        
     
     create_patch_mesh( four_quads, 9, coords, 4, indices, QUADRILATERAL, err );
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
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     pd.set_mesh_set(mesh_set1);


     double coords[] = { 1,.5, 0,
                         0, 0, 0,
                         1, 0, 0,
                         2, 0, 0,
                         2, 1, 0,
                         2, 2, 0,
                         1, 2, 0,
                         0, 2, 0,
                         0, 1, 0,
                         3, 0, 0,
                         3, 1, 0,
                         3, 2, 0 };

     size_t indices[] = { 1,  2,  0, 8,
                          2,  3,  4, 0,
                          8,  0,  6, 7,
                          0,  4,  5, 6,
                          3,  9, 10, 4,
                          4, 10, 11, 5 };
     
     create_patch_mesh( pd, 12, coords, 6, indices, QUADRILATERAL, err );

     MsqVertex* vert_array = pd.get_vertex_array(err);
     vert_array[1].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[2].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[3].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[5].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[6].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[7].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[8].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[9].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[10].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
     vert_array[11].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED);
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
     double coords[] = { 1, 1, -1,
                         0, 0, -1,
                         1, 0, -1,
                         2, 0, -1,
                         2, 1, -1,
                         2, 2, -1,
                         1, 2, -1,
                         0, 2, -1,
                         0, 1, -1,
                         3, 0, -1,
                         3, 1, -1,
                         3, 2, -1,

                         1,.5, 0,
                         0, 0, 0,
                         1, 0, 0,
                         2, 0, 0,
                         2, 1, 0,
                         2, 2, 0,
                         1, 2, 0,
                         0, 2, 0,
                         0, 1, 0,
                         3, 0, 0,
                         3, 1, 0,
                         3, 2, 0,

                         1, 1, 1,
                         0, 0, 1,
                         1, 0, 1,
                         2, 0, 1,
                         2, 1, 1,
                         2, 2, 1,
                         1, 2, 1,
                         0, 2, 1,
                         0, 1, 1,
                         3, 0, 1,
                         3, 1, 1,
                         3, 2, 1 };
     
     size_t connectivity[] = { 1, 2, 0, 8, 13, 14, 12, 20, // 0
                               2, 3, 4, 0, 14, 15, 16, 12, // 1
                               8, 0, 6, 7, 20, 12, 18, 19, // 2
                               0, 4, 5, 6, 12, 16, 17, 18, // 3
                               3, 9, 10, 4, 15, 21, 22, 16, // 4
                               4, 10, 11, 5, 16, 22, 23, 17, // 5
                               13, 14, 12, 20, 25, 26, 24, 32, // 6
                               14, 15, 16, 12, 26, 27, 28, 24, // 7
                               20, 12, 18, 19, 32, 24, 30, 31, // 8
                               12, 16, 17, 18, 24, 28, 29, 30, // 9
                               15, 21, 22, 16, 27, 33, 34, 28, // 10
                               16, 22, 23, 17, 28, 34, 35, 29 }; // 11

     create_patch_mesh( pd, 36, coords, 12, connectivity, HEXAHEDRON, err );
     
     MsqVertex* verts = pd.get_vertex_array(err); MSQ_CHKERR(err);

     verts[0].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[1].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[2].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[3].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[4].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[5].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[6].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[7].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[8].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[9].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[10].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[11].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
      
     verts[12].set_vertex_flag(MsqVertex::MSQ_NO_VTX_FLAG); 
     verts[13].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[14].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[15].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[16].set_vertex_flag(MsqVertex::MSQ_NO_VTX_FLAG); 
     verts[17].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[18].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[19].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[20].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[21].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[22].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[23].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
      
     verts[24].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[25].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[26].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[27].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[28].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[29].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[30].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[31].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[32].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[33].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[34].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
     verts[35].set_vertex_flag(MsqVertex::MSQ_HARD_FIXED); 
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
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     triPatch.set_mesh_set(mesh_set1);

     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         0.5, sqrt(3.0)/2.0, 0.0,
                         2.0, -4.0, 2.0 };
     

     const size_t conn[] = { 0, 1, 2, 0, 3, 1 };
     
     create_patch_mesh( triPatch, 4, coords, 2, conn, TRIANGLE, err );
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
     PlanarDomain* msq_geom = new PlanarDomain(s_norm, pnt);
     mesh_set1->set_domain_constraint(msq_geom, err); MSQ_CHKERR(err);
     quadPatch.set_mesh_set(mesh_set1);

     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         1.0, 1.0, 0.0,
                         0.0, 1.0, 0.0,
                         2.0, -1.0, .5,
                         1.5, 1.0, 1.0 };
     
     const size_t conn[] = { 0, 1, 2, 3, 1, 4, 5, 2 };
     
     create_patch_mesh( quadPatch, 6, coords, 2, conn, QUADRILATERAL, err );
   }
  
     /* Patch used in several quality metric tests.
        Our tet patch is made of two tets.  tet_1 is a perfect
        equilateral (the ideal for most metrics).  tet_2 is an arbitrary
        tet.
     */
   inline void create_qm_two_tet_patch(PatchData &tetPatch, MsqError &err)
   {
     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         0.5, sqrt(3.0)/2.0, 0.0,
                         0.5, sqrt(3.0)/6.0, sqrt(2.0)/sqrt(3.0),
                         2.0, 3.0, -.5 };
     

     const size_t conn[] = { 0, 1, 2, 3, 1, 4, 2, 3 };
     
     create_patch_mesh( tetPatch, 5, coords, 2, conn, TETRAHEDRON, err );
   }
   /* Patch used in seveal quality metric tests.
      Our hex patch is made of two hexes.  hex_1 is a perfect
      unit cube (the ideal for most metrics).  hex_2 is an arbitrary
      hex.
   */
   inline void create_qm_two_hex_patch(PatchData &hexPatch, MsqError &err)
   {
     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         1.0, 1.0, 0.0,
                         0.0, 1.0, 0.0,
                         0.0, 0.0, 1.0,
                         1.0, 0.0, 1.0,
                         1.0, 1.0, 1.0,
                         0.0, 1.0, 1.0,
                         2.0, 0.0, 0.0,
                         2.0, 1.0, 0.0,
                         2.0,-1.0, 1.0,
                         3.0, 2.0, 1.0 };
     
     const size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7,
                             1, 8, 9, 2, 5, 10, 11, 6 };
                             
     create_patch_mesh( hexPatch, 12, coords, 2, conn, HEXAHEDRON, err ); 
   }
   
} // namespace

#endif // PatchDataInstances_hpp
