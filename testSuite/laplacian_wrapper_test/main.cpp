// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 14-Nov-02 at 17:33:19 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifdef USE_STD_INCLUDES
#include <iostream>
#else
#include <iostream.h>
#endif

#ifdef USE_C_PREFIX_INCLUDES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include "Mesquite.hpp"
#include "TSTT_Base.h"
#include "MesquiteUtilities.hpp" //  for writeVtkMesh()
#include "MesquiteError.hpp"
#include "MeshSet.hpp"

// algorythms
#include "LaplacianIQ.hpp"


using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "main"
int main()
{     
  char file_name[128];
  /* Reads a TSTT Mesh file */
  TSTT::Mesh_Handle mesh;
  TSTT::MeshError tstt_err;
  TSTT::Mesh_Create(&mesh, &tstt_err);
  strcpy(file_name, "../../meshFiles/2D/VTK/square_quad_2.vtk");
  TSTT::Mesh_Load(mesh, file_name, &tstt_err);
  
  // Mesquite error object
  MsqError err;
  
  // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
  // creates an intruction queue
  LaplacianIQ laplacian_smoother;

  

  writeVtkMesh("original_mesh", mesh, err); MSQ_CHKERR(err);
  
  // launches optimization on mesh_set1
  laplacian_smoother.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
  writeVtkMesh("smoothed_mesh", mesh, err); MSQ_CHKERR(err);

}
