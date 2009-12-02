#include "ShapeImprovementWrapper.hpp"
#include "MeshImpl.hpp"
#include "CLArgs.hpp"
#include "MsqError.hpp"

#include "domain.hpp"

#include <iostream>

using namespace Mesquite;

int main( int argc, char* argv[] )
{
  const double zero = 0.0;
  double default_cpu_time = 0.0;
  CLArgs::DoubleRangeArg l2_norm( &zero, 0 ), cpu_time( default_cpu_time, &zero, 0 );
  
  CLArgs args( "msqshape",
               "Run Shape Improvement smoother for input mesh.",
               "Read VTK file, smooth, and re-write file." );
  args.double_flag( 'n', "GradL2Norm", "termination graident L2 norm", &l2_norm );
  args.double_flag( 't', "Cpu Seconds", "time-out", &cpu_time );
                    
  add_domain_args( args );
  args.add_required_arg( "input_file" );
  args.add_required_arg( "output_file" );

  std::vector<std::string> files;
  if (!args.parse_options( argc, argv, files, std::cerr )) {
    args.print_usage( std::cerr );
    exit(1);
  }
  std::string input_file = files[0];
  std::string output_file = files[1];
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( input_file.c_str(), err );
  if (err) {
    std::cerr << "ERROR READING FILE: " << input_file << std::endl
                    << err << std::endl;
    return 2;
  }
  MeshDomain* domain = process_domain_args( &mesh );

  if (l2_norm.seen()) {
    ShapeImprovementWrapper smoother( err, cpu_time.value(), l2_norm.value() );
    if (err) {
      std::cerr << "Error constructing smoother" << std::endl
                      << err << std::endl;
      return 2;
    }
    smoother.run_instructions( &mesh, domain, err );
  }
  else {
    ShapeImprovementWrapper smoother( err, cpu_time.value() );
    if (err) {
      std::cerr << "Error constructing smoother" << std::endl
                      << err << std::endl;
      return 2;
    }
    smoother.run_instructions( &mesh, domain, err );
  }
  if (err) {
    std::cerr << err << std::endl;
    return 3;
  }
  
  mesh.write_vtk( output_file.c_str(), err );
  if (err) {
    std::cerr << "ERROR WRITING FILE: " << output_file << std::endl
                    << err << std::endl;
    return 2;
  }
  
  return 0;
}

