// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 23-May-03 at 18:04:38 by Thomas Leurent
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
#include "MeshImpl.hpp"
#include "MesquiteError.hpp"

using namespace Mesquite;

using std::cout;
using std::endl;

#undef __FUNC__
#define __FUNC__ "main"
int main(int argc, char* argv[])
{
  Mesquite::MsqError err;
  char in_file_name[256];
  char out_file_name[256];
  double OF_value = 1.;
  
  // command line arguments
  if (argc!=3){
    cout << "Input meshfile name needed as first argument.\n"
      "Output meshfile name needed as second argument.\n" << endl;
    return -1;
  }
  else{
    cout << " given 2 command line arguments.\n";
    strcpy(in_file_name, argv[1]);
    strcpy(out_file_name, argv[2]);
  }
#ifndef MSQ_USING_EXODUS
  err.set_msg("Exodus not enabled in this build of Mesquite");
  MSQ_CHKERR(err);
  return -1;
#else
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  cout<<"\nReading VTK file.\n";
  mesh->read_vtk(in_file_name, err);MSQ_CHKERR(err);
  cout<<"Writing Exodus file.\n";
  mesh->write_exodus(out_file_name,err);MSQ_CHKERR(err);
#endif
}
