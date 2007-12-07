/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file main.cpp
 *  \brief Test syncronous boundary case
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "LPtoPTemplate.hpp"
#include "InstructionQueue.hpp"
#include "ConjugateGradient.hpp"
#include "FeasibleNewton.hpp"
#include "MeshImpl.hpp"
#include "QualityAssessor.hpp"

#include "data.hpp"

#include <iostream>
#include <cstdlib>

using namespace Mesquite;
const char default_out_file[] = "synchronous.vtk";
double default_x = 0.25;

void usage( const char* argv0, bool brief = true )
{
  std::ostream& str = brief ? std::cerr : std::cout;
  
  str << "Usage: " << argv0 
      << " [-x <coord>]"
      << " [-j|-n|-d]"
      << " [-r|-c]"
      << " [<output_file>]"
      << std::endl;
  if (brief) {
    str << "       " << argv0 << " -h" << std::endl;
    std::exit(1);
  }
  
  str << "  -x  Specify X coordinate value for mesh (default is " << default_x << ")" << std::endl
      << "  -j  Use ConjugateGradient solver (default)" << std::endl
      << "  -n  Use FeasibleNewton solver" << std::endl
      << "  -r  Use IdealWeightInverseMeanRation metric" << std::endl
      << "  -c  Use ConditionNumber metric (default)" << std::endl
      << "Default output file is \"" << default_out_file << '"' << std::endl;
  
  std::exit(0);
}

char mSolver = '\0', mMetric = '\0';
const char* outputFile = default_out_file;
double input_x = default_x;

void parse_options( char* argv[], int argc )
{
  bool next_arg_is_x = false;
  for (int i = 1; i < argc; ++i) {
    if (next_arg_is_x) {
      next_arg_is_x = false;
      char* endptr = 0;
      input_x = std::strtod( argv[i], &endptr );
      if (endptr && *endptr)
        usage(argv[0]);
      continue;
    }
      
    if (argv[i][0] != '-') {
      if (outputFile != default_out_file)
        usage(argv[0]);
      outputFile = argv[i];
      continue;
    }
    
    for (const char* p = argv[i]+1; *p; ++p) {
      switch (*p) {
        case 'x': 
          next_arg_is_x = true; 
          break;
        
        case 'j':
        case 'n':
          if (mSolver) 
            usage(argv[0]);
          mSolver = *p;
          break;
        
        case 'r':
        case 'c':
          if (mMetric) 
            usage(argv[0]);
          mMetric = *p;
          break;
 
        default:
          usage( argv[0], *p != 'h' );
          break;
      }
    }
  }
    
  if (next_arg_is_x)
    usage(argv[0]);
  
    // default values
  if (!mMetric)
    mMetric = 'c';
  if (!mSolver)
    mSolver = 'j';
}

int main( int argc, char* argv[] )
{
  parse_options( argv, argc );
  
  MeshImpl mesh;
  MyDomain domain;
  MsqError err;
  
  create_input_mesh( input_x, mesh, err );
  if (MSQ_CHKERR(err)) { std::cerr << err << std::endl; return 2; }
  
  domain.setup( &mesh, err );
  if (MSQ_CHKERR(err)) { std::cerr << err << std::endl; return 2; }
  
  QualityMetric* metric = 0;
  if (mMetric == 'c')
    metric = new ConditionNumberQualityMetric;
  else
    metric = new IdealWeightInverseMeanRatio;
  
  LPtoPTemplate function( 1, metric );
  
  VertexMover* solver = 0;
  if (mSolver == 'j')
    solver = new ConjugateGradient( &function );
  else
    solver = new FeasibleNewton( &function );
    
  if (PatchSetUser* psu = dynamic_cast<PatchSetUser*>(solver)) 
    psu->use_global_patch();
  
  TerminationCriterion inner;
  inner.add_criterion_type_with_double( TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, 1e-4, err );
  inner.write_vtk_timesteps( "synchronous" );
  solver->set_inner_termination_criterion( &inner );  
  
  InstructionQueue q;
  QualityAssessor qa( metric, 10 );
  q.add_quality_assessor( &qa, err );
  q.set_master_quality_improver( solver, err );
  q.add_quality_assessor( &qa, err );
  q.run_instructions( &mesh, &domain, 0, err );
  delete solver;
  delete metric;
  
  if (MSQ_CHKERR(err)) 
    { std::cerr << err << std::endl; return 3; }
    
  mesh.write_vtk( outputFile, err );
  if (MSQ_CHKERR(err)) 
    { std::cerr << err << std::endl; return 2; }
  
  return 0;
}

  
