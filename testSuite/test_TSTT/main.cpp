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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 10-Feb-04 at 22:44:58 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include "Mesquite.hpp"
#include "MeshTSTT.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "MeshWriter.hpp"

// algorithms
#include "IdealWeightInverseMeanRatio.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "SmartLaplacianSmoother.hpp"

#ifndef MSQ_USE_OLD_IO_HEADERS
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#else
#include <iostream.h>
#endif

#ifdef MSQ_USE_OLD_C_HEADERS
#include <cstdlib>
#else
#include <stdlib.h>
#endif


#if   defined( MSQ_TSTT_USE_MOAB )     && MSQ_TSTT_USE_MOAB
# include "TSTTM_MOAB_MoabMesh.hh"
  typedef TSTTM_MOAB::MoabMesh ImplType;
#elif defined( MSQ_TSTT_USE_OVERTURE ) && MSQ_TSTT_USE_OVERTURE
# include "TSTT_Overture_Mesh.hh"
  typedef TSTT_Overture::Mesh ImplType;
#elif defined( MSQ_TSTT_USE_AOMD )     && MSQ_TSTT_USE_AOMD
# include "TSTT_LocalTSTTMesh_Impl.hh"
  typedef TSTT::LocalTSTTMesh ImplType;
#else
# error
#endif

#include "TSTTB.hh"
#include "TSTTM.hh"


using namespace Mesquite;

const char* const default_file_name = "../../meshFiles/3D/VTK/large_box_hex_1000.vtk";

void usage()
{
  cout << "main [-N] [filename]" << endl;
  cout << "  -N : Use native representation instead of TSTT implementation\n";
  cout << "  If no file name is specified, will use \"" 
       << default_file_name << '"' << endl;
  exit (1);
}

  // Construct a MeshTSTT from the file
Mesh* get_tstt_mesh( const char* file_name );

  // Construct a MeshImpl from the file
Mesh* get_native_mesh( const char* file_name );

  // Run FeasibleNewton solver
int run_global_smoother( Mesh* mesh, MsqError& err );

  // Run SmoothLaplacian solver
int run_local_smoother( Mesh* mesh, MsqError& err );

int main(int argc, char* argv[])
{
  Mesquite::MsqPrintError err(cout);
  
    // command line arguments
  const char* file_name = 0;
  bool use_native = false, opts_done = false;
  for (int arg = 1; arg < argc; ++arg)
  {
    if (!opts_done && argv[arg][0] == '-')
    {
      if (!strcmp( argv[arg], "-N"))
        use_native = true;
      else if(!strcmp( argv[arg], "--"))
        opts_done = true;
      else
        usage();
    }
    else if (!file_name)
      file_name = argv[arg];
    else
      usage();
  }
  if (!file_name)
  {
    file_name = default_file_name;
    cout << "No file specified: using default: " << default_file_name << endl;
  }  
  
    // Try running a global smoother on the mesh
  Mesh* mesh = use_native ? 
               get_native_mesh(file_name) : 
               get_tstt_mesh(file_name);
  MeshWriter::write_vtk(mesh, "original.vtk", err); 
  if (err) return 1;
  cout << "Wrote \"original.vtk\"" << endl;
  run_global_smoother( mesh, err );
  if (err) return 1;
  
    // Try running a local smoother on the mesh
  mesh = use_native ? 
         get_native_mesh(file_name) : 
         get_tstt_mesh(file_name);
  run_local_smoother( mesh, err );
  if (err) return 1;
  
  return 0;
}


int run_global_smoother( Mesh* mesh, MsqError& err )
{
  double OF_value = 0.0001;

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::SUM, err); 
  if (err) return 1;
  
  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  
  // creates the feas newt optimization procedures
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  if (err) return 1;
  
  QualityAssessor stop_qa=QualityAssessor(mean_ratio,QualityAssessor::AVERAGE, err);
  if (err) return 1;
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_double(
           TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, OF_value, err);
  if (err) return 1;
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err) return 1;
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;
   
  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;
  
  queue1.add_quality_assessor(&stop_qa, err); 
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(mesh, err);
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "feasible-newton-result.vtk", err); 
  if (err) return 1;
  cout << "Wrote \"feasible-newton-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}

int run_local_smoother( Mesh* mesh, MsqError& err )
{
  double OF_value = 0.0001;

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  ShapeQualityMetric* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::SUM, err); 
  if (err) return 1;
  
  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  
  // creates the smart laplacian optimization procedures
  SmartLaplacianSmoother* pass1 = new SmartLaplacianSmoother( obj_func, err );
  if (err) return 1;
  
  QualityAssessor stop_qa=QualityAssessor(mean_ratio,QualityAssessor::AVERAGE,err);
  if (err) return 1;
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_double(
           TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, OF_value, err);
  if (err) return 1;
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err) return 1;
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;
   
  // adds 1 pass of pass1 to mesh_set
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;
  
  queue1.add_quality_assessor(&stop_qa, err); 
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(mesh, err);
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "smart-laplacian-result.vtk", err); 
  if (err) return 1;
  cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}
 


Mesh* get_tstt_mesh( const char* file_name )
{
  
  ImplType mesh_instance = ImplType::_create();
  TSTTM::Mesh tstt_mesh = mesh_instance;

  tstt_mesh.load(0, file_name);
  
  // **************************************************************
  // THIS IS A HACK -- READ FIXED FLAG FROM INPUT FILE DIRECTLY
  // AND SET TAG IN TSTT MESH IMPLEMENTATION
  //  - assumes TSTTM implementation does not already set tag
  //  - assumes order of vertex handles in TSTTM implementation 
  //    is the same as the vertex order in the file. THIS IS A
  //    BIG ASSUMPTION!  It works for the current MOAB 
  //    implementation.
  // **************************************************************
  FILE* file = fopen( file_name, "r" );
  char buffer[256];
  int count = 0;
  while (fgets( buffer, sizeof(buffer), file ))
    if (sscanf(buffer, "POINT_DATA %d", &count) == 1)
      break;
    // Found fixed flags?
  if (count)
  {
      // skip two lines
    fgets( buffer, sizeof(buffer), file );
    fgets( buffer, sizeof(buffer), file );
      // read flags from file
    std::vector<int> flags(count);
    for (std::vector<int>::iterator iter = flags.begin(); iter != flags.end(); ++iter)
      if (fscanf( file, "%d", &*iter ) != 1) {
        fprintf(stderr, "Error reading fixed flag from \"%s\"\n", file_name );
        fclose( file );
        exit (2);
      }
    fclose( file );
      // Get tag interface
    TSTTB::ArrTag tag_iface = tstt_mesh;
    if (!tag_iface) {
      fprintf(stderr, "TSTTM impelementation does not provide tag array interface\n");
      exit (2);
    }
      // Store flags in tag on vertices
    try {
        // Get vertices
      int num_vtx;
      sidl::array<void*> vertices;
      tstt_mesh.getEntities( tstt_mesh.getRootSet(),
                             TSTTM::EntityType_VERTEX,
                             TSTTM::EntityTopology_POINT,
                             vertices, num_vtx );
      if (num_vtx != count) {
        fprintf(stderr, "File \"%s\" contains fixed flags for %d vertices, "
                        "while TSTTM implementation has %d vertices.\n",
                        file_name, count, num_vtx );
        exit(2);
      }
        // Create tag
      TagHandle tag_handle;
      tag_iface.createTag( VERTEX_FIXED_TAG_NAME, 
                           1, 
                           TSTTB::TagValueType_INTEGER, 
                           tag_handle );
        // Set tag data on vertices
      sidl::array<int> values;
      int32_t lower = 0, upper = flags.size(), stride = 1;
      values.borrow( &flags[0], 1, &lower, &upper, &stride );
      tag_iface.setIntArrData( vertices, num_vtx, tag_handle, values, flags.size() );
    } 
    catch( ... ) {
      std::cerr << "Error setting fixed tag on vertices\n";
      exit(2);
    } 
  }
  // **************************************************************
  // END UGLY HACK
  // **************************************************************
  

  MsqError err;
  Mesquite::MeshTSTT* mesh = MeshTSTT::create( tstt_mesh, tstt_mesh.getRootSet(), err );
  if (err)
  {
    cerr << err << endl;
    exit (1);
  }
  
  return mesh;
}
  


Mesh* get_native_mesh( const char* file_name )
{
  MsqError err;
  MeshImpl* mesh = new MeshImpl;
  mesh->read_vtk( file_name, err );
  if (err)
  {
    cerr << err << endl;
    exit(3);
  }
  
  return mesh;
}


