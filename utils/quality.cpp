#include "CLArgs.hpp"
#include "QualityAssessor.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "SizeMetric.hpp"
#include "TMPQualityMetric.hpp"
#include "IdealTargetCalculator.hpp"
#include "Target3DShape.hpp"
#include "Target2DShape.hpp"
#include "InstructionQueue.hpp"
#include "MsqError.hpp"
#include "MeshImpl.hpp"
#include "QuadLagrangeShape.hpp"
#include "TetLagrangeShape.hpp"
#include "TriLagrangeShape.hpp"
#include "ElementMaxQM.hpp"

using namespace Mesquite;

#include <vector>
#include <algorithm>

void tag_fixed_elements( MeshImpl& mesh, const char* tag_name );

int main( int argc, char* argv[] )
{
  int two;
  CLArgs::ToggleArg freeonly;
  CLArgs::IntRangeArg histogram( &two );
  CLArgs args( "msqquality", "Assess mesh quality",
               "Caculate various quality metrics for a mesh," 
               "and optinally export a VTK file for which quality "
               "values are stored as attribute data." );
  args.int_flag( 'H', "ints", "Print histograms with specified number of intervals", &histogram );
  args.toggle_flag( 'f', "Assess quality only for elements with at least one free vertex", &freeonly );
  args.add_required_arg( "input_file" );
  args.add_optional_arg( "output_file" );
  msq_std::vector<msq_std::string> files;
  if (!args.parse_options( argc, argv, files, msq_stdio::cerr )) {
    args.print_usage( msq_stdio::cerr );
    return 1;
  }
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( files.front().c_str(), err );
  if (err) {
    msq_stdio::cerr << err << msq_stdio::endl 
                    << "Failed to read file: " << files.front() << msq_std::endl;
    return 2;
  }

  QualityAssessor qa( true, freeonly.value(), "INVERTED" );
  IdealWeightInverseMeanRatio imr;
  SizeMetric size;
  IdealTargetCalculator tc;
  Target3DShape tm_3d;
  Target2DShape tm_2d;
  TMPQualityMetric tmp( &tc, &tm_2d, &tm_3d );
  ElementMaxQM max_tmp( &tmp );
  
  int intervals = histogram.seen() ? histogram.value() : 0;
  qa.add_quality_assessment( &imr, intervals, 0.0, "InverseMeanRatio" );
  qa.add_quality_assessment( &size, intervals, 0.0, "Size" );
  qa.add_quality_assessment( &max_tmp, intervals, 0.0, "TMP_Shape" );

  QuadLagrangeShape quad;
  TriLagrangeShape tri;
  TetLagrangeShape tet;
  InstructionQueue q;
  q.set_mapping_function( &quad );
  q.set_mapping_function( &tri );
  q.set_mapping_function( &tet );
  
  q.add_quality_assessor( &qa, err );
  q.run_instructions( &mesh, err );
  if (err) {
    msq_stdio::cerr << err << msq_stdio::endl;
    return 3;
  }
  
  if (files.size() > 1) {
    tag_fixed_elements( mesh, "FIXED_ELEMS" );
    mesh.write_vtk( files[1].c_str(), err );
    if (err) {
      msq_stdio::cerr << err << msq_stdio::endl 
                      << "Failed to write file: " << files[1] << msq_std::endl;
      return 2;
    }
  }
  
  return 0;
}

#define TFEERR(A) if (A) { \
     msq_stdio::cerr << "Failed to tag fixed elements" << msq_stdio::endl; \
     msq_stdio::cerr << err << msq_stdio::endl; \
     return; }

void tag_fixed_elements( MeshImpl& mesh, const char* tag_name )
{
  MsqError err;
  
  msq_std::vector<Mesh::ElementHandle> elements;
  mesh.get_all_elements( elements, err );
  TFEERR(err);
  msq_std::vector<int> tag_vals(elements.size());
  
  msq_std::vector<Mesh::VertexHandle> verts;
  msq_std::vector<size_t> junk;
  bool fixed[64];
  int count = 0;
  for (size_t i = 0; i < elements.size(); ++i) {
    verts.clear(); junk.clear();
    mesh.elements_get_attached_vertices( &elements[i], 1, verts, junk, err );
    TFEERR(err);
    mesh.vertices_get_fixed_flag( &verts[0], fixed, verts.size(), err );
    TFEERR(err);
    
    bool* end = fixed+verts.size();
    if (end == std::find(fixed, end, false)) {
      ++count;
      tag_vals[i] = 1;
    }
    else 
      tag_vals[i] = 0;
  }
  
  if (!count)
    return;
    
  TagHandle tag = mesh.tag_create( tag_name, Mesh::INT, 1, 0, err );
  TFEERR(err);
  mesh.tag_set_element_data( tag, elements.size(), &elements[0], &tag_vals[0], err );
  TFEERR(err);  
  msq_stdio::cerr << "Counted " << count << " fixed elements" << msq_stdio::endl;
}

