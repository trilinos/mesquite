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
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 13-Apr-04 at 10:57:52
//  LAST-MOD: 18-Oct-04 by J.Kraftcheck
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include "meshfiles.h"

#ifndef MSQ_USE_OLD_IO_HEADERS
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#endif

#ifdef MSQ_USE_OLD_C_HEADERS
#include <cstdlib>
#else
#include <stdlib.h>
#endif


#include "Mesquite.hpp"
#include "MeshImpl.hpp"
//#include "MeshMunson.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "MeshWriter.hpp"
#include "TargetReader.hpp"
#include "TargetWriter.hpp"
#include "ReferenceMesh.hpp"
#include "UnitWeight.hpp"

// algorithms
#include "I_DFT.hpp"
#include "ConcreteTargetCalculators.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
using namespace Mesquite;

void usage()
{
  cerr << "Usage: main [filename] [objective function val]" << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  MsqPrintError err(cout);

  const char* file_name = MESH_FILES_DIR "2D/VTK/tfi_horse10x4-14.vtk";
  double OF_value = 0.;
  
  // command line arguments
  if (argc == 1)
  {
    cerr << "Warning: No file specified, using default: " << file_name << endl;
  }
  
  if (argc > 1)
  {
    file_name = argv[1];
  }
  if (argc > 2)
  {
    char* end_ptr;
    OF_value = strtod( argv[2], &end_ptr );
    if (!*argv[2] || *end_ptr)
      usage();
  }
  if (argc > 3)
  {
    usage();
  }
  
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(file_name, err);
  if (err) return 1;
  
  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  TargetReader reader(true);
  UnitWeight weights;
  I_DFT mean_ratio( &reader, &weights );

  Mesquite::MeshImpl ref_mesh;
  ref_mesh.read_vtk(MESH_FILES_DIR "2D/VTK/tfi_horse10x4-12.vtk", err);
  if (err) return 1;
  Mesquite::ReferenceMesh rm(&ref_mesh);
  DeformingDomainGuides843 target( &rm );
  TargetWriter writer( DistanceFromTarget::get_dft_sample_pts(), &target );
  queue1.add_target_calculator( &writer, err ); 
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate obj_func(&mean_ratio, 1, err);
  if (err) return 1;
  
  // creates the steepest descentfeas newt optimization procedures
//  ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  FeasibleNewton pass1( &obj_func, true );
  pass1.use_global_patch();
  
  QualityAssessor stop_qa( &mean_ratio );
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  if (OF_value!=0) {
    tc_inner.add_criterion_type_with_double(
           TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, OF_value, err);
    if (err) return 1;
    pass1.set_inner_termination_criterion(&tc_inner);
  }
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  if (err) return 1;
  pass1.set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err); 
  if (err) return 1;
   
  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(&pass1, err); 
  if (err) return 1;
  
  queue1.add_quality_assessor(&stop_qa, err); 
  if (err) return 1;

  MeshWriter::write_gnuplot( &ref_mesh, "ref_mesh", err );
  if (err) return 1;

  MeshWriter::write_gnuplot( &mesh, "ori_mesh", err );
  if (err) return 1;
  
  // launches optimization on mesh_set1
  queue1.run_instructions( &mesh, err ); 
  if (err) return 1;
  
  MeshWriter::write_gnuplot( &mesh, "smo_mesh", err );
  if (err) return 1;

  print_timing_diagnostics( cout );
  return 0;
}
