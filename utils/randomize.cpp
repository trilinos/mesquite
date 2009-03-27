#include "Randomize.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"
#include "MeshImpl.hpp"

#include "MeshDomain1D.hpp"
#include "PlanarDomain.hpp"
#include "CylinderDomain.hpp"
#include "SphericalDomain.hpp"
#include "DomainClassifier.hpp"

#include "PatchData.hpp"
#include "MsqVertex.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "PMeanPTemplate.hpp"
#include "TerminationCriterion.hpp"
#include <assert.h>

#include "CLArgs.hpp"

#include <iostream>
#include <iomanip>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

using namespace Mesquite;


const char SPHERE_FLAG = 'S';
const char PLANE_FLAG = 'P';
const char CYLINDER_FLAG = 'C';
const char LINE_FLAG = 'l';
const char CIRCLE_FLAG = 'c';
const char POINT_FLAG = 'v';
const char INVALID_FLAG = 'i';
const char PERCENT_FLAG = 'p';
const char SKIN_FLAG = 's';
const char UNOPTIMIZE_FLAG = 'u';

class SphereDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  SphereDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool SphereDomainArg::value( const msq_std::vector<double>& list )
{
  double rad = list[0];
  if (rad <= 0.0)
    return false;
  Vector3D center(0,0,0);
  if (list.size() == 4)
    center.set( list[1], list[2], list[3] );
  domList.push_back( new SphericalDomain( center, rad ) );
  dimList.push_back( 2 );
  return true;
}

class CylinderDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  CylinderDomainArg( msq_std::vector<MeshDomain*>& domlist,
                     msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool CylinderDomainArg::value( const msq_std::vector<double>& vals )
{
  double rad = vals[0];
  Vector3D normal( vals[1], vals[2], vals[3] );
  Vector3D point(0,0,0);
  if (vals.size() == 7)
    point.set( vals[4], vals[5], vals[6] );
  domList.push_back( new CylinderDomain( rad, normal, point ) );
  dimList.push_back( 2 );
  return true;
}

class PlanarDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  PlanarDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool PlanarDomainArg::value( const msq_std::vector<double>& list )
{
  Vector3D normal( list[0], list[1], list[2] );
  Vector3D point(0,0,0);
  if (list.size() == 6)
    point.set( list[3], list[4], list[5] );
  domList.push_back( new PlanarDomain( normal, point ) );
  dimList.push_back( 2 );
  return true;
}

class LineDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  LineDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool LineDomainArg::value( const msq_std::vector<double>& vals )
{
  Vector3D dir( vals[0], vals[1], vals[2] );
  Vector3D point(0,0,0);
  if (vals.size() == 6)
    point.set( vals[3], vals[4], vals[5] );
  LineDomain* pdom = new LineDomain( point, dir );
  domList.push_back( pdom );
  dimList.push_back( 1 );
  return true;
}

class CircleDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  CircleDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool CircleDomainArg::value( const msq_std::vector<double>& vals )
{
  double rad = vals[0];
  Vector3D normal( vals[1], vals[2], vals[3] );
  Vector3D point(0,0,0);
  if (vals.size() == 7)
    point.set( vals[4], vals[5], vals[6] );
  CircleDomain* pdom = new CircleDomain( point, normal, rad );
  domList.push_back( pdom );
  dimList.push_back( 1 );
  return true;
}

class PointDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  PointDomainArg( msq_std::vector<MeshDomain*>& domlist,
                  msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool PointDomainArg::value( const msq_std::vector<double>& vals )
{
  Vector3D point( vals[0], vals[1], vals[2] );
  PointDomain* pdom = new PointDomain( point );
  domList.push_back( pdom );
  dimList.push_back( 0 );
  return true;
}


class UnOptimizer : public VertexMover
{
public:
  
  UnOptimizer( ObjectiveFunction* of ) : objectiveFunction(of) {}
  
  virtual ~UnOptimizer() {}
  
  virtual msq_std::string get_name() const;
  
  virtual PatchSet* get_patch_set();
  
protected:

    virtual void initialize(PatchData &pd, MsqError &err);
    
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
                                         
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    
    virtual void cleanup();

private:
  
    msq_std::vector<size_t> adjVtxList;
    VertexPatches patchSet;
    ObjectiveFunction* objectiveFunction;
};

msq_std::string UnOptimizer::get_name() const { return "UnOptimize"; }
PatchSet* UnOptimizer::get_patch_set() { return &patchSet; }
void UnOptimizer::initialize( PatchData&, MsqError& ) {}
void UnOptimizer::initialize_mesh_iteration(PatchData& , MsqError& ) {}
void UnOptimizer::terminate_mesh_iteration(PatchData& , MsqError& ) {}
void UnOptimizer::cleanup() {}

void UnOptimizer::optimize_vertex_positions( PatchData &pd, MsqError &err) {
  assert( pd.num_free_vertices() == 1 && pd.vertex_by_index(0).is_free_vertex() );
  msq_std::vector<Vector3D> grad(1);
  double val, junk, coeff;
  bool state;
  
  state = objectiveFunction->evaluate_with_gradient( ObjectiveFunction::CALCULATE,
                                                     pd, val, grad, err );
  MSQ_ERRRTN(err);
  if (!state) {
    MSQ_SETERR(err)(MsqError::INVALID_MESH);
    return;
  }
  grad[0] /= grad[0].length();
  
  PatchDataVerticesMemento* memento = pd.create_vertices_memento( err ); MSQ_ERRRTN(err);
  msq_std::auto_ptr<PatchDataVerticesMemento> deleter( memento );
  pd.get_minmax_edge_length( junk, coeff );
  
  for (int i = 0; i < 100; ++i) {
    pd.set_free_vertices_constrained( memento, &grad[0], 1, coeff, err ); MSQ_ERRRTN(err);
    state = objectiveFunction->evaluate( ObjectiveFunction::CALCULATE, pd, val, true, err );
    MSQ_ERRRTN(err);
    if (state)
      break;
    coeff *= 0.5;
  }
  if (!state) {
    pd.set_to_vertices_memento( memento, err );
  }
}


int main( int argc, char* argv[] )
{
  const double default_fraction = 0.05;
  const double zero = 0.0;
  int one = 1;
  CLArgs::ToggleArg allow_invalid( false );
  CLArgs::ToggleArg skin_mesh( false );
  CLArgs::DoubleRangeArg rand_percent( default_fraction, &zero, 0 );
  CLArgs::IntRangeArg unoptimize( 0, &one, 0 );
  
  msq_std::vector<MeshDomain*> domains;
  msq_std::vector<int> domain_dims;
  SphereDomainArg     sphere_arg( domains, domain_dims );
  CylinderDomainArg cylinder_arg( domains, domain_dims );
  PlanarDomainArg      plane_arg( domains, domain_dims );
  CircleDomainArg     circle_arg( domains, domain_dims );
  LineDomainArg         line_arg( domains, domain_dims );
  PointDomainArg       point_arg( domains, domain_dims );
  
  const char* SPHERE_VALUES[] = { "rad", "x", "y", "z" };
  const char* CYLINDER_VALUES[] = { "rad", "i", "j", "k", "x", "y", "z" };
  CLArgs args( "vtkrandom",
               "Randomize mesh vertex locations.",
               "Read VTK file, randomize locations of containded vertices, and re-write file." );
  args.toggle_flag( INVALID_FLAG, "Allow inverted elements in output", &allow_invalid );
  args.toggle_flag( SKIN_FLAG, "Mark boundary vertices as fixed (default if no domain specified)", &skin_mesh );
  args.double_flag( PERCENT_FLAG, "fract", "Randomize fraction", &rand_percent );
  args.int_flag( UNOPTIMIZE_FLAG, "N", "Use UnOptimizer with N passes rather than Randomize", &unoptimize );
  args.double_list_flag( SPHERE_FLAG, "Spherical domain as center and radius", &sphere_arg );
  args.limit_list_flag( SPHERE_FLAG, 4, SPHERE_VALUES );
  args.limit_list_flag( SPHERE_FLAG, 1, SPHERE_VALUES );
  args.double_list_flag( PLANE_FLAG, "Planar domain as normal and point", &plane_arg );
  args.limit_list_flag( PLANE_FLAG, 3, CYLINDER_VALUES+1 );
  args.limit_list_flag( PLANE_FLAG, 6, CYLINDER_VALUES+1 );
  args.double_list_flag( CYLINDER_FLAG, "Cylindrical domain as point, axis, and radius", &cylinder_arg );
  args.limit_list_flag( CYLINDER_FLAG, 4, CYLINDER_VALUES );
  args.limit_list_flag( CYLINDER_FLAG, 7, CYLINDER_VALUES );
  args.double_list_flag( LINE_FLAG, "Linear domain as direction and point", &line_arg );
  args.limit_list_flag( LINE_FLAG, 3, CYLINDER_VALUES+1 );
  args.limit_list_flag( LINE_FLAG, 6, CYLINDER_VALUES+1 );
  args.double_list_flag( CIRCLE_FLAG, "Circular domain as radius, normal, and center", &circle_arg );
  args.limit_list_flag( CIRCLE_FLAG, 4, CYLINDER_VALUES );
  args.limit_list_flag( CIRCLE_FLAG, 7, CYLINDER_VALUES );
  args.double_list_flag( POINT_FLAG, "Point domain", &point_arg );
  args.limit_list_flag( POINT_FLAG, 3, SPHERE_VALUES+1 );
  args.add_required_arg( "input_file" );
  args.add_required_arg( "output_file" );

  msq_std::vector<msq_std::string> files;
  if (!args.parse_options( argc, argv, files, msq_stdio::cerr )) {
    args.print_usage( msq_stdio::cerr );
    exit(1);
  }
  msq_std::string input_file = files[0];
  msq_std::string output_file = files[1];
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( input_file.c_str(), err );
  if (err) {
    msq_stdio::cerr << "ERROR READING FILE: " << input_file << msq_stdio::endl
                    << err << msq_stdio::endl;
    return 2;
  }
  
  MeshDomain* domain = 0;
  DomainClassifier combined_domain;
  if (!domains.empty()) {
    domain = &combined_domain;
    DomainClassifier::classify_geometrically( combined_domain,
                                              &mesh,
                                              1e-4,
                                              &domains[0],
                                              &domain_dims[0],
                                              domains.size(),
                                              err );
    if (err) {
      msq_stdio::cerr << err << msq_stdio::endl;
      return 3;
    }
  }
  
  if (domains.empty() || skin_mesh.value() ) {
    mesh.mark_skin_fixed( err, false );
    if (err) {
      msq_stdio::cerr << err << msq_stdio::endl;
      return 3;
    }
  }    
  
  TerminationCriterion tc;
  QualityAssessor qa( false );
  InstructionQueue q;
  Randomize op( rand_percent.value() );
  IdealWeightInverseMeanRatio metric;
  PMeanPTemplate of( 1, &metric );
  UnOptimizer op2( &of );
  if (unoptimize.seen()) {
    tc.add_iteration_limit( unoptimize.value() );
    op2.set_outer_termination_criterion( &tc );
    q.add_preconditioner( &op, err );
    q.set_master_quality_improver( &op2, err );
  }
  else {
    q.set_master_quality_improver( &op, err );
  }
  q.add_quality_assessor( &qa, err );
  q.run_instructions( &mesh, domain, err );
  if (err) {
    msq_stdio::cerr << err << msq_stdio::endl;
    return 3;
  }

  int inverted, junk;
  if (qa.get_inverted_element_count( inverted, junk, err ) && inverted ) {
    if (allow_invalid.value())
      msq_stdio::cerr << "Warning: output mesh contains " << inverted << " inverted elements" << msq_stdio::endl;
    else {
      msq_stdio::cerr << "Error: output mesh contains " << inverted << " inverted elements" << msq_stdio::endl;
      return 4;
    }
  }
  
  mesh.write_vtk( output_file.c_str(), err );
  if (err) {
    msq_stdio::cerr << "ERROR WRITING FILE: " << output_file << msq_stdio::endl
                    << err << msq_stdio::endl;
    return 2;
  }
  
  return 0;
}

