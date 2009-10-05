/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Target2DSurfOrientTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "DomainSurfaceOrientation.hpp"
#include "TargetSurfaceOrientation.hpp"
#include "TMPQualityMetric.hpp"
#include "InstructionQueue.hpp"
#include "PMeanPTemplate.hpp"
#include "CompositeOFScalarMultiply.hpp"
#include "CompositeOFAdd.hpp"
#include "ConjugateGradient.hpp"
#include "MeshImpl.hpp"
#include "PlanarDomain.hpp"
#include "SphericalDomain.hpp"
#include "MeshWriter.hpp"
#include "MsqVertex.hpp"
#include "Settings.hpp"
#include "IdealTargetCalculator.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "ReferenceMesh.hpp"
#include "QualityAssessor.hpp"
#include "TerminationCriterion.hpp"

#include "Target2DShapeSize.hpp"
#include "Target2DShapeSizeOrient.hpp"
typedef Mesquite::Target2DShapeSize Metric2D_Type_Ideal;
typedef Mesquite::Target2DShapeSizeOrient Metric2D_Type_Ref;


#include "UnitUtil.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "meshfiles.h"
#include <iostream>
#include <vector>

using namespace Mesquite;
using namespace std;

const bool DISABLE_SPHERE_OUTPUT_FILES = false;
const bool DISABLE_PLANAR_OUTPUT_FILES = true;
const double METRIC_3D_ORIENT_FACTOR = 2.0;
const double METRIC_3D_ORIENT_P = 1.0;
const double METRIC_2D_PLANE_P = 1.0;

class Target2DSurfOrientTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(Target2DSurfOrientTest);

  CPPUNIT_TEST (compare_none_plane_xy_to_yz);
  CPPUNIT_TEST (compare_target_plane_xy_to_yz);
  CPPUNIT_TEST (compare_domain_plane_xy_to_yz);

  CPPUNIT_TEST (compare_none_plane_xy_to_xz);
  CPPUNIT_TEST (compare_target_plane_xy_to_xz);
  CPPUNIT_TEST (compare_domain_plane_xy_to_xz);

  CPPUNIT_TEST (compare_none_plane_xy_to_neg_xy);
  CPPUNIT_TEST (compare_target_plane_xy_to_neg_xy);
  CPPUNIT_TEST (compare_domain_plane_xy_to_neg_xy);
  
  CPPUNIT_TEST (compare_none_target_plane_xy);
  CPPUNIT_TEST (compare_target_domain_plane_xy);

  CPPUNIT_TEST (try_none_ideal_sphere);
  CPPUNIT_TEST (try_none_ref_sphere);

  CPPUNIT_TEST (compare_target_domain_ideal_sphere);
  CPPUNIT_TEST (compare_target_domain_ref_sphere);

  CPPUNIT_TEST_SUITE_END();

  MeshImpl myMesh, refMesh;
  std::vector<MsqVertex> savedCoords;
  
  void load_plane_mesh();
  void load_sphere_mesh();
  void rotate_mesh( PlanarDomain::Plane from,
                    PlanarDomain::Plane to );

  /**\param metric_3D 0->no_3d_orientation, 1->TargetSurfaceOrientation, 2->DomainSurfaceOrientation
   * \param ref_mesh true->RefMeshTargetCalculator, false->IdealTargetCalculator
   */
  void smooth( int metric_3D, bool ref_mesh, MeshDomain* dom );
  
  void store_coordinates();
  void compare_coordinates( double epsilon );

  void compare_rotation( PlanarDomain::Plane target_plane,
                         int metric_3D, const char* output_file = 0 );
  void compare_neg_xy( int metric_3D, const char* output_file = 0 );
                         
  void compare_metrics( int metric1, int metric2, const char* output = 0 );

  void smooth_sphere( int metric, bool ref_mesh, const char* output = 0 );

public:

  void compare_none_plane_xy_to_yz() 
    { compare_rotation( PlanarDomain::YZ, 0, "none_xy_to_yz" ); }
  void compare_target_plane_xy_to_yz()
    { compare_rotation( PlanarDomain::YZ, 1, "target_xy_to_yz" ); }
  void compare_domain_plane_xy_to_yz()
    { compare_rotation( PlanarDomain::YZ, 2, "domain_xy_to_yz" ); }

  void compare_none_plane_xy_to_xz()
    { compare_rotation( PlanarDomain::XZ, 0, "none_xy_to_xz" ); }
  void compare_target_plane_xy_to_xz()
    { compare_rotation( PlanarDomain::XZ, 1, "target_xy_to_xz" ); }
  void compare_domain_plane_xy_to_xz()
    { compare_rotation( PlanarDomain::XZ, 2, "domain_xy_to_xz" ); }

  void compare_none_plane_xy_to_neg_xy() 
    { compare_neg_xy( 0, "none_xy_to_neg_xy" ); }
  void compare_target_plane_xy_to_neg_xy()
    { compare_neg_xy( 1, "target_xy_to_neg_xy" ); }
  void compare_domain_plane_xy_to_neg_xy()
    { compare_neg_xy( 2, "domain_xy_to_neg_xy" ); }
  
  void try_none_ideal_sphere() 
    { smooth_sphere( 0, false, "none_ideal" ); }
  
  void try_none_ref_sphere() 
    { smooth_sphere( 0, true, "none_ref" ); }
  
  void compare_none_target_plane_xy()
    { compare_metrics( 0, 1, "none_vs_target" ); }
  void compare_target_domain_plane_xy()
    { compare_metrics( 1, 2, "target_vs_domain" ); }

  void compare_target_domain_ideal_sphere() 
  {
    smooth_sphere( 1, false, "target_ideal" );
    store_coordinates();
    smooth_sphere( 2, false, "domain_ideal" );
    compare_coordinates( 1e-3 );
  }
  
  void compare_target_domain_ref_sphere()
  {
    smooth_sphere( 1, true, "target_ref" );
    store_coordinates();
    smooth_sphere( 2, true, "domain_ref" );
      // result meshes are expected to be roughly the same,
      // but not exactly, so use large tolerance when comparing
    compare_coordinates( 0.2 );
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Target2DSurfOrientTest, "Target2DSurfOrientTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Target2DSurfOrientTest, "Regression");

void Target2DSurfOrientTest::compare_rotation( PlanarDomain::Plane target_plane, 
                                               int metric_3D, 
                                               const char* output_file )
{
  if (output_file)
    cout << endl << "****************************" << output_file << "**************************" << endl << endl;
  MsqPrintError err(cerr);
  PlanarDomain xy(PlanarDomain::XY);
  PlanarDomain pl(target_plane);

  load_plane_mesh();
  smooth( metric_3D, false, &xy );
  if (!DISABLE_PLANAR_OUTPUT_FILES && output_file) {
    string n( output_file );
    n += ".expected";
    MeshWriter::write_gnuplot( &myMesh, n.c_str(), err );
    ASSERT_NO_ERROR(err);
  }
  store_coordinates();
  
  load_plane_mesh();
  rotate_mesh( PlanarDomain::XY, target_plane );
  smooth( metric_3D, false, &pl );
  rotate_mesh( target_plane, PlanarDomain::XY );
  if (!DISABLE_PLANAR_OUTPUT_FILES && output_file) {
    string n( output_file );
    n += ".actual";
    MeshWriter::write_gnuplot( &myMesh, n.c_str(), err );
    ASSERT_NO_ERROR(err);
  }
  compare_coordinates( 1e-3 );
}

void Target2DSurfOrientTest::compare_neg_xy( int metric_3D, const char* output_file )
{
  if (output_file)
    cout << endl << "****************************" << output_file << "**************************" << endl << endl;
  MsqPrintError err(cerr);
  PlanarDomain xy(PlanarDomain::XY);
  PlanarDomain pl(Vector3D(0,0,-1),Vector3D(0,0,0));
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<MsqVertex> coords;
  size_t i;

  load_plane_mesh();
  smooth( metric_3D, false, &xy );
  if (!DISABLE_PLANAR_OUTPUT_FILES && output_file) {
    string n( output_file );
    n += ".expected";
    MeshWriter::write_gnuplot( &myMesh, n.c_str(), err );
    ASSERT_NO_ERROR(err);
  }
  store_coordinates();
  
  load_plane_mesh();

    // get coordinates
  myMesh.get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!vertices.empty());
  coords.resize( vertices.size() );
  myMesh.vertices_get_coordinates( &vertices[0],  &coords[0], vertices.size(), err );
  ASSERT_NO_ERROR(err);
  
    // rotate and set coordinates
  for (i = 0; i < vertices.size(); ++i) {
    coords[i][0] = -coords[i][0];
    myMesh.vertex_set_coordinates( vertices[i], coords[i], err );
    ASSERT_NO_ERROR(err);
  }
  
  smooth( metric_3D, false, &pl );
  
    // rotate back for comparison
  myMesh.vertices_get_coordinates( &vertices[0],  &coords[0], vertices.size(), err );
  ASSERT_NO_ERROR(err);
  for (i = 0; i < vertices.size(); ++i) {
    coords[i][0] = -coords[i][0];
    myMesh.vertex_set_coordinates( vertices[i], coords[i], err );
    ASSERT_NO_ERROR(err);
  }
  
  if (!DISABLE_PLANAR_OUTPUT_FILES && output_file) {
    string n( output_file );
    n += ".actual";
    MeshWriter::write_gnuplot( &myMesh, n.c_str(), err );
    ASSERT_NO_ERROR(err);
  }
  compare_coordinates( 1e-3 );
}

void Target2DSurfOrientTest::compare_metrics( int metric1, int metric2, const char* output )
{
  if (output)
    cout << endl << "****************************" << output << "**************************" << endl << endl;
  MsqPrintError err(cerr);
  PlanarDomain dom( PlanarDomain::XY );

  load_plane_mesh();
  smooth( metric1, false, &dom );
  if (!DISABLE_PLANAR_OUTPUT_FILES && output) {
    string n( output );
    n += ".expected";
    MeshWriter::write_gnuplot( &myMesh, n.c_str(), err );
    ASSERT_NO_ERROR(err);
  }
  store_coordinates();
  
  load_plane_mesh();
  smooth( metric2, false, &dom );
  if (!DISABLE_PLANAR_OUTPUT_FILES && output) {
    string n( output );
    n += ".actual";
    MeshWriter::write_gnuplot( &myMesh, n.c_str(), err );
    ASSERT_NO_ERROR(err);
  }
  compare_coordinates( 1e-3 );
}


void Target2DSurfOrientTest::smooth_sphere( int metric, bool ref_mesh, const char* output )
{
  if (output)
    cout << endl << "****************************" << output << "**************************" << endl << endl;
  MsqPrintError err(cerr);
  SphericalDomain dom( Vector3D(0.0, 0.0, 0.0), 5.0 );
  load_sphere_mesh();
  smooth( metric, ref_mesh, &dom );
  if (!DISABLE_SPHERE_OUTPUT_FILES && output) {
    string n(output);
    n += ".vtk";
    myMesh.write_vtk( n.c_str(), err );
  }
}

void Target2DSurfOrientTest::load_plane_mesh()
{
  MsqPrintError err(cerr);
  size_t i;
  
    // load the mesh 
  myMesh.read_vtk( MESH_FILES_DIR "2D/VTK/square_quad_10_rand.vtk", err );
  ASSERT_NO_ERROR(err);
  
    // move mesh into XY plane
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<MsqVertex> coords;
  myMesh.get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!vertices.empty());
  coords.resize( vertices.size() );
  myMesh.vertices_get_coordinates( &vertices[0],  &coords[0], vertices.size(), err );
  ASSERT_NO_ERROR(err);
  for (i = 0; i < vertices.size(); ++i) {
    coords[i][2] = 0.0;
    myMesh.vertex_set_coordinates( vertices[i], coords[i], err );
    ASSERT_NO_ERROR(err);
  }
}

void Target2DSurfOrientTest::load_sphere_mesh()
{
  MsqPrintError err(cerr);
  myMesh.read_vtk( MESH_FILES_DIR "3D/VTK/quad_sphere.vtk", err );
  ASSERT_NO_ERROR(err);
  refMesh.read_vtk( MESH_FILES_DIR "3D/VTK/quad_sphere_ref.vtk", err );
  ASSERT_NO_ERROR(err);
}

void Target2DSurfOrientTest::rotate_mesh( PlanarDomain::Plane from,
                                          PlanarDomain::Plane to )
{
  MsqPrintError err(cerr);
  size_t i;
  
    // get coordinates
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<MsqVertex> coords;
  myMesh.get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!vertices.empty());
  coords.resize( vertices.size() );
  myMesh.vertices_get_coordinates( &vertices[0],  &coords[0], vertices.size(), err );
  
    // rotate and set coordinates
  for (i = 0; i < vertices.size(); ++i) {
    switch (from) {
      case PlanarDomain::XY:
        if (to == PlanarDomain::YZ) {
          coords[i][2] = -coords[i][0];
          coords[i][0] = 0.0;
        }
        else if (to == PlanarDomain::XZ) {
          coords[i][2] = -coords[i][1];
          coords[i][1] = 0.0;
        }
        break;
      case PlanarDomain::YZ:
        if (to == PlanarDomain::XY) {
          coords[i][0] = -coords[i][2];
          coords[i][2] = 0.0;
        }
        else if (to == PlanarDomain::XZ) {
          coords[i][0] = coords[i][1];
          coords[i][1] = 0.0;
        }
        break;
      case PlanarDomain::XZ:
        if (to == PlanarDomain::XY) {
          coords[i][1] = -coords[i][2];
          coords[i][2] = 0.0;
        }
        else if (to == PlanarDomain::YZ) {
          coords[i][1] = coords[i][0];
          coords[i][0] = 0.0;
        }
        break;
    }
    myMesh.vertex_set_coordinates( vertices[i], coords[i], err );
    ASSERT_NO_ERROR(err);
  }
}

void Target2DSurfOrientTest::smooth( int metric_3D, bool ref_mesh, MeshDomain* dom )
{
  MsqPrintError err(cerr);
  
  IdealTargetCalculator ideal_target( !ref_mesh );
  ReferenceMesh rm( &refMesh );
  RefMeshTargetCalculator ref_target( &rm );

  TargetCalculator* tc = ref_mesh ? (TargetCalculator*)&ref_target : (TargetCalculator*)&ideal_target;

  Metric2D_Type_Ideal ideal_tmetric;
  Metric2D_Type_Ref ref_tmetric;
  TargetMetric2D* tmetric = ref_mesh ? (TargetMetric2D*)&ref_tmetric : (TargetMetric2D*)&ideal_tmetric;
  TMPQualityMetric qmetric(tc, tmetric, NULL);
  TargetSurfaceOrientation ometric1(tc);
  DomainSurfaceOrientation ometric2;
  
  QualityMetric* ometric = metric_3D == 2 ? (QualityMetric*)&ometric2 : (QualityMetric*)&ometric1;
  PMeanPTemplate pmean_3D( METRIC_3D_ORIENT_P, ometric );
  CompositeOFScalarMultiply scale_3D( METRIC_3D_ORIENT_FACTOR, &pmean_3D );
  ObjectiveFunction* OF_3D = METRIC_3D_ORIENT_FACTOR == 1.0 ? (ObjectiveFunction*)&pmean_3D : (ObjectiveFunction*)&scale_3D;

  PMeanPTemplate pmean_2D( METRIC_2D_PLANE_P, &qmetric );
  CompositeOFAdd sum_OF( OF_3D, &pmean_2D );
  ObjectiveFunction* OF = metric_3D ? (ObjectiveFunction*)&sum_OF : (ObjectiveFunction*)&pmean_2D;
  
  TerminationCriterion inner;
  inner.add_absolute_vertex_movement( 1e-3 );
  inner.add_iteration_limit( 20 );
  
  ConjugateGradient solver( OF );
  solver.use_global_patch();
  solver.set_inner_termination_criterion( &inner );
  
  QualityAssessor assessor;
  assessor.add_quality_assessment( &qmetric );
  assessor.add_quality_assessment( ometric );
  
  InstructionQueue q;
  q.add_quality_assessor( &assessor, err );
  q.set_master_quality_improver( &solver, err );
  q.add_quality_assessor( &assessor, err );
  
  q.run_instructions( &myMesh, dom, err );
  ASSERT_NO_ERROR(err);
}

void Target2DSurfOrientTest::store_coordinates()
{
  MsqPrintError err(cerr);
  std::vector<Mesh::VertexHandle> vertices;
  myMesh.get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!vertices.empty());
  savedCoords.resize( vertices.size() );
  myMesh.vertices_get_coordinates( &vertices[0],  &savedCoords[0], vertices.size(), err );
  ASSERT_NO_ERROR(err);
}

void Target2DSurfOrientTest::compare_coordinates( double epsilon )
{
  MsqPrintError err(cerr);
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<MsqVertex> coords;
  myMesh.get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!vertices.empty());
  coords.resize( vertices.size() );
  myMesh.vertices_get_coordinates( &vertices[0],  &coords[0], vertices.size(), err );
  ASSERT_NO_ERROR(err);
  
  CPPUNIT_ASSERT_EQUAL( coords.size(), savedCoords.size() );
  for (size_t i = 0; i < coords.size(); ++i)
    CPPUNIT_ASSERT_VECTORS_EQUAL( savedCoords[i], coords[i], epsilon );
}

