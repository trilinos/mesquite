// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 21-Aug-02 at 16:04:29
//  LAST-MOD: 31-Oct-02 at 13:05:35 by Thomas Leurent
//
// RCS Infos:
// ==========
//    Author: $Author$
//        Id: $Id$
//  Revision: $Revision$
//      Date: $Date$
//    locker: $Locker$
//
//
// DESCRIPTION:
// ============
/*! \file MesquiteUtilities.cpp

describe MesquiteUtilities.cpp here

\author Thomas Leurent -- tleurent@mcs.anl.gov
*/
// DESCRIP-END.
//

#include <stdio.h>
#include "MesquiteUtilities.hpp"
#include "MsqMeshEntity.hpp"
#include <assert.h>
//Michael do not check this in
//#include <iostream.h>
using namespace Mesquite; 


#undef __FUNC__
#define __FUNC__ "writeVtkMesh"
/*! \fn writeVtkMesh(const char filebase[128], TSTT::cMesh_Handle mesh_h, MsqError &err)
  
    This function takes a TSTT mesh handle and writes a mesh in VTK format to filebase.vtk.
    The VTK format can deal with hybrid and non-hybrid meshes of triangles, quads,
    tetrahedrons and hexahedrons. 
*/
void Mesquite::writeVtkMesh(const char filebase[128], TSTT::cMesh_Handle mesh_h,
                                   MsqError &err)
{
  int num_vertices=0;
  int num_elem = 0;
  char vtk_file[128];
  FILE *vtk_fp;
  TSTT::MeshError tstt_err = 0;
  int* id_1 = new int;
  int* id_2 = new int;

  sprintf(vtk_file,"%s.%s",filebase,"vtk");

  if ((vtk_fp = fopen(vtk_file,"w")) == NULL) {
    err.set_msg("cannot open vtk_file."); MSQ_CHKERR(err);
    return;
  }

  // retrieves all the vertices
  TSTT::Entity_Handle *vtx;
  TSTT::Mesh_GetEntities(mesh_h, TSTT::VERTEX, &vtx, &num_vertices, &tstt_err);
  assert(!tstt_err);

  // retrieves all the elements. First it tries to retrieve the TSTT::REGIONs (Tet/Hex)
  // and if there isn't any, it retrieves the TSTT::FACEs (Tri/Quad).
  TSTT::Entity_Handle *elements;
  TSTT::Mesh_GetEntities(mesh_h, TSTT::REGION, &elements, &num_elem, &tstt_err);
  if (num_elem == 0) {
    TSTT::Mesh_GetEntities(mesh_h, TSTT::FACE, &elements, &num_elem, &tstt_err);
    assert(!tstt_err);
  }
  // if we still don't have any element, sets an error.
  if (num_elem == 0) {
    err.set_msg("no elements found in the TSTT::Mesh_Handle."); MSQ_CHKERR(err);
  }
  
  
  fprintf(vtk_fp,"# vtk DataFile Version 2.0\n");
  fprintf(vtk_fp,"Mesquite Mesh %s .\n", filebase);
  fprintf(vtk_fp,"ASCII\n");
  fprintf(vtk_fp,"DATASET UNSTRUCTURED_GRID\n");


  // gets and prints out the coordinates of each vertex
  fprintf(vtk_fp,"POINTS %d float\n", num_vertices);
  int num_coords = 3*num_vertices;
  double* vtx_coords = new double[num_coords];
  TSTT::Entity_GetVertexCoords(mesh_h,
                               (TSTT::cEntity_Handle *) vtx,
                               num_vertices, TSTT::INTERLEAVED,
                               &vtx_coords, &num_coords,
                               &tstt_err);
  assert(!tstt_err);
  for (int i=0; i<num_vertices; ++i) {
    fprintf(vtk_fp,"%f   %f   %f\n", vtx_coords[3*i], vtx_coords[3*i+1], vtx_coords[3*i+2]);
  }


  //checks the topology of the elements
  TSTT::EntityTopology *topologies = new TSTT::EntityTopology[num_elem];
  TSTT::Entity_GetTopology(mesh_h, (TSTT::cEntity_Handle *) elements, 
                     num_elem, topologies, &tstt_err);
  assert(!tstt_err);
  int cells_table_size =0; // this is an info needed by VTK. see p.600.
  bool includes_tris = false;
  bool includes_quads = false;
  bool includes_tets = false;
  bool includes_hexs = false;
  for (int i=0; i<num_elem; ++i) {
    cells_table_size += 1 + MsqMeshEntity::vertex_count(topologies[i], err); MSQ_CHKERR(err);
    if (topologies[i] == TSTT::TRIANGLE)       includes_tris=true;
    if (topologies[i] == TSTT::QUADRILATERAL)  includes_quads=true;
    if (topologies[i] == TSTT::TETRAHEDRON)    includes_tets=true;
    if (topologies[i] == TSTT::HEXAHEDRON)     includes_hexs=true;
  }
  // We only allow tri/quad and tet/hex meshes ...
  if ( includes_tris && (includes_tets || includes_tets) ||
       includes_quads && (includes_tets || includes_tets) ) {
    err.set_msg("invalid mix of elements in the mesh. Only tri/quad and tet/hex allowed.");
    return;
  }
  fprintf(vtk_fp,"CELLS %d %d\n", num_elem, cells_table_size);
  
  // get and prints out the connectivity for each element
  TSTT::Entity_Handle* elem_vertices = new TSTT::Entity_Handle[8];
  for (int e=0; e<num_elem; ++e) {
    int* csr_pointer; // dummy
    int* csr_data; // dummy
    int num_elem_vtx; num_elem_vtx=8;
    TSTT::Entity_GetAdjacencies( mesh_h, (TSTT::cEntity_Handle *) &(elements[e]),
                                 1, TSTT::VERTEX, &elem_vertices,
                                 &csr_pointer, &csr_data, 
                                 &num_elem_vtx, &tstt_err );
    assert(!tstt_err);
    int index=0;
    int num_vtx_in_element = MsqMeshEntity::vertex_count(topologies[e], err);
    fprintf(vtk_fp,"%d ", num_vtx_in_element); MSQ_CHKERR(err);
    // for each vertex in the element
    for (int i=0; i<num_vtx_in_element; ++i)
    {
        // finds out the vertex index in the local scope
        // list given its global TSTT ID
      bool proper_index = false;
      index=0;
      // gets global ID of triangle vertex
      TSTT::Entity_GetGlobalID(mesh_h, (TSTT::cEntity_Handle *) &(elem_vertices[i]),
                               1, &id_1, &tstt_err);
      while (!proper_index)
      {
          // gets global ID of vertex at position "index" in our local scope list
        TSTT::Entity_GetGlobalID(mesh_h, (TSTT::cEntity_Handle *) &(vtx[index]),
                                 1, &id_2, &tstt_err);
          // compare those IDs
        if (*id_1 == *id_2)
          proper_index = true;
        if (index>=num_vertices)
        {
          err.set_msg("index pb"); MSQ_CHKERR(err); return;
        }
        ++index;
      }
        // print the index 
      fprintf(vtk_fp,"%d ", index-1); // 0-base for VTK
    }
    fprintf(vtk_fp,"\n");
  }

  // writes the topology of each element.
  fprintf(vtk_fp, "CELL_TYPES %d\n", num_elem);
  for (int e=0; e<num_elem; ++e) {
    // converts from TSTT to VTK topology type.
    int vtk_type;
    switch (topologies[e]) {
    case TSTT::TRIANGLE:
      vtk_type = 5;
      break;
    case TSTT::QUADRILATERAL:
      vtk_type = 9;
      break;
    case TSTT::TETRAHEDRON:
      vtk_type = 10;
      break;
    case TSTT::HEXAHEDRON:
      vtk_type = 12;
      break;
    default:
      break;
    }
    fprintf(vtk_fp, "%d\n", vtk_type);
  }
  
  // writes which points are fixed
  fprintf(vtk_fp, "POINT_DATA %d\n", num_vertices);
  fprintf(vtk_fp, "SCALARS fixed float\n");
  fprintf(vtk_fp, "LOOKUP_TABLE default\n");

  // retrieves value of tag "fixed" .
  char bnd_tag_name[6];
  strcpy(bnd_tag_name, "fixed");
  void* bnd_tag_handle;
  TSTT::Mesh_tagGetHandle (mesh_h, bnd_tag_name, &bnd_tag_handle, &tstt_err);
  assert(!tstt_err);
  for (int i=0; i<num_vertices; ++i) {
    int* on_boundary;
    int tag_size;
    TSTT::Mesh_GetTag_Entity(mesh_h,
                             (TSTT::cEntity_Handle) vtx[i],
                             bnd_tag_handle, (void**)&on_boundary, &tag_size, &tstt_err);
    assert(!tstt_err);
    fprintf(vtk_fp, "%d\n", *on_boundary);
  }
  

  delete id_1;
  delete id_2;
  delete[] vtx_coords;
  delete[] topologies;
  delete[] elem_vertices;
}


#undef __FUNC__
#define __FUNC__ "writeTSTTShowMeMesh"
/* This is mostly a hack ... not foolproof ... be sure to hand it a triangular mesh */
void Mesquite::writeTSTTShowMeMesh(const char filebase[128], TSTT::cMesh_Handle mesh_h,
                                   MsqError &err)
{
  int num_vertices=0;
  int num_tri = 0;
  char node_file[128], tri_file[128];
  FILE *node_fp, *tri_fp;
  TSTT::MeshError tstt_err = 0;
  int* id_1 = new int;
  int* id_2 = new int;

  sprintf(node_file,"%s.%s",filebase,"node");
  sprintf(tri_file,"%s.%s",filebase,"ele");

  if ((node_fp =  fopen(node_file,"w")) == NULL) {
    printf("cannot open node_file\n");    exit(0);
  }
  if ((tri_fp =  fopen(tri_file,"w")) == NULL) {
    printf("cannot open tri_file\n");    exit(0);
  }

  TSTT::Entity_Handle *vtx;
  TSTT::Mesh_GetEntities(mesh_h, TSTT::VERTEX, &vtx, &num_vertices, &tstt_err);
  assert(!tstt_err);
  TSTT::Entity_Handle *tri;
  TSTT::Mesh_GetEntities(mesh_h, TSTT::FACE, &tri, &num_tri, &tstt_err);
  assert(!tstt_err);
    
  fprintf(node_fp,"%d 2 0 1\n", num_vertices);
  fprintf(tri_fp,"%d %d %d\n", num_tri,3,0);

  // gets and prints out the coordinates
  int num_coords = 3*num_vertices;
  double* vtx_coords = new double[num_coords];
  TSTT::Entity_GetVertexCoords(mesh_h,
                               (TSTT::cEntity_Handle *) vtx,
                               num_vertices, TSTT::INTERLEAVED,
                               &vtx_coords, &num_coords,
                               &tstt_err);
  assert(!tstt_err);
  for (int i=0; i<num_vertices; ++i) {
    fprintf(node_fp,"%d   %f   %f   %f\n",i+1, vtx_coords[3*i],
            vtx_coords[3*i+1], vtx_coords[3*i+2]);
  }

  // get and prints out the connectivity for each triangle
  TSTT::Entity_Handle* tri_vertices = new TSTT::Entity_Handle[3];
  for (int e=0; e<num_tri; ++e) {
    int* csr_pointer; // dummy
    int* csr_data; // dummy
    int num_tri_vtx=3;
    TSTT::Entity_GetAdjacencies( mesh_h, (TSTT::cEntity_Handle *) &(tri[e]),
                                 1, TSTT::VERTEX, &tri_vertices,
                                 &csr_pointer, &csr_data, 
                                 &num_tri_vtx, &tstt_err );
    assert(!tstt_err);
    int index=0;
    fprintf(tri_fp,"%d   ",e+1); // prints triangle number.
    // for each vertex in the triangle
    for (int i=0; i<3; ++i)
    {
        // finds out the vertex index in the local scope
        // list given its global TSTT ID
      bool proper_index = false;
      index=0;
      // gets global ID of triangle vertex
      TSTT::Entity_GetGlobalID(mesh_h, (TSTT::cEntity_Handle *) &(tri_vertices[i]),
                               1, &id_1, &tstt_err);
      while (!proper_index)
      {
          // gets global ID of vertex at position "index" in our local scope list
        TSTT::Entity_GetGlobalID(mesh_h, (TSTT::cEntity_Handle *) &(vtx[index]),
                                 1, &id_2, &tstt_err);
          // compare those IDs
        if (*id_1 == *id_2)
          proper_index = true;
        if (index>=num_vertices)
        {
          err.set_msg("index pb"); MSQ_CHKERR(err); return;
        }
        ++index;
      }
        // print the index 
      fprintf(tri_fp,"%d   ", index);
    }
    fprintf(tri_fp,"\n");
  }
  delete id_1;
  delete id_2;
  delete[] tri_vertices;
}


#undef __FUNC__
#define __FUNC__ "writeTSTTFacetMesh"
/* This is mostly a hack ... not foolproof ... be sure to hand it a triangular mesh */
void Mesquite::writeTSTTFacetMesh(const char filebase[128], TSTT::cMesh_Handle mesh_h,
                                   MsqError &err)
{
  int num_vertices=0;
  int num_tri = 0;
  char facet_file[128];
  FILE *facet_fp;
  TSTT::MeshError tstt_err = 0;
  int* id_1 = new int;
  int* id_2 = new int;

    //sprintf(node_file,"%s.%s",filebase,"node");
  sprintf(facet_file,"%s.%s",filebase,"facet");

  if ((facet_fp =  fopen(facet_file,"w")) == NULL) {
    printf("cannot open facet_file\n");    exit(0);
  }

  TSTT::Entity_Handle *vtx;
  TSTT::Mesh_GetEntities(mesh_h, TSTT::VERTEX, &vtx, &num_vertices, &tstt_err);
  assert(!tstt_err);
  TSTT::Entity_Handle *tri;
  TSTT::Mesh_GetEntities(mesh_h, TSTT::FACE, &tri, &num_tri, &tstt_err);
  assert(!tstt_err);
    /*
      fprintf(facet_fp,"%d %d\n", num_nodes, num_tri);
  for (i=0;i<num_nodes;i++) {
    fprintf(facet_fp,"%d   %15.10f   %15.10f   %15.10f\n",i+1,mesh->vtx[i]->coord[XDIR],
  	    mesh->vtx[i]->coord[YDIR], mesh->vtx[i]->coord[ZDIR]);
  }
  for (i=0;i<num_tri;i++) {
    fprintf(facet_fp,"%d   %d   %d   %d\n",i+1,mesh->tri[i]->vtx_id[0]+1,
            mesh->tri[i]->vtx_id[1]+1,mesh->tri[i]->vtx_id[2]+1);
  }
    */
  fprintf(facet_fp,"%d %d\n", num_vertices, num_tri);

  // gets and prints out the coordinates
  int num_coords = 3*num_vertices;
  double* vtx_coords = new double[num_coords];
  TSTT::Entity_GetVertexCoords(mesh_h,
                               (TSTT::cEntity_Handle *) vtx,
                               num_vertices, TSTT::INTERLEAVED,
                               &vtx_coords, &num_coords,
                               &tstt_err);
  assert(!tstt_err);
  for (int i=0; i<num_vertices; ++i) {
    fprintf(facet_fp,"%d   %f   %f   %f\n",i+1, vtx_coords[3*i],
            vtx_coords[3*i+1], vtx_coords[3*i+2]);
  }

  // get and prints out the connectivity for each triangle
  TSTT::Entity_Handle* tri_vertices = new TSTT::Entity_Handle[3];
  for (int e=0; e<num_tri; ++e) {
    int* csr_pointer; // dummy
    int* csr_data; // dummy
    int num_tri_vtx=3;
    TSTT::Entity_GetAdjacencies( mesh_h, (TSTT::cEntity_Handle *) &(tri[e]),
                                 1, TSTT::VERTEX, &tri_vertices,
                                 &csr_pointer, &csr_data, 
                                 &num_tri_vtx, &tstt_err );
    assert(!tstt_err);
    int index=0;
    fprintf(facet_fp,"%d   ",e+1); // prints triangle number.
    // for each vertex in the triangle
    for (int i=0; i<3; ++i) {
      // finds out the vertex index in the local scope list given its global TSTT ID
      bool proper_index = false;
      index=0;
      // gets global ID of triangle vertex
      TSTT::Entity_GetGlobalID(mesh_h, (TSTT::cEntity_Handle *) &(tri_vertices[i]),
                               1, &id_1, &tstt_err);
      while (!proper_index) {
        // gets global ID of vertex at position "index" in our local scope list
        TSTT::Entity_GetGlobalID(mesh_h, (TSTT::cEntity_Handle *) &(vtx[index]),
                                 1, &id_2, &tstt_err);
        // compare those IDs
        if (*id_1 == *id_2)
          proper_index = true;
        if (index>=num_vertices) {
          err.set_msg("index pb"); MSQ_CHKERR(err); return; }
        ++index;
      }
      // print the index 
      fprintf(facet_fp,"%d   ", index);
    }
    fprintf(facet_fp,"\n", index);
  }
  delete id_1;
  delete id_2;
  delete[] tri_vertices;
}


