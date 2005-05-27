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
    kraftche@cae.wisc.edu   
   
  ***************************************************************** */

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <iostream.h>
#else
#  include <iostream>
   using std::cout;
   using std::endl;
#endif

#include <stdlib.h>


#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "LPtoPTemplate.hpp"

// algorithms
#include "IdealWeightMeanRatio.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "ConjugateGradient.hpp"
using namespace Mesquite;

// Use CPPUNIT_ASSERT in code so it's easy to convert to a unit test later.
#define CPPUNIT_ASSERT(A) \
  do { if (!(A)) { \
  msq_stdio::cout << "Assertion Failed: " << #A << msq_stdio::endl; \
  msq_stdio::cout << "  File: " << __FILE__ << msq_stdio::endl; \
  msq_stdio::cout << "  Line: " << __LINE__ << msq_stdio::endl; \
  return true; \
  } } while (false)

// Given a mesh with a single free vertex located at the origin,
// move the vertex to the specified position, smooth the mesh,
// and verify that the vertex was moved back to the origin by
// the smoother.
bool smooth_mesh( Mesh* mesh,
                  Mesh::VertexHandle free_vertex_at_origin, 
                  Vector3D initial_free_vertex_position,
                  QualityMetric* metric );

int main()
{
  Mesquite::MsqPrintError err(cout);
  QualityMetric* metrics[] = { new IdealWeightMeanRatio,
                               new IdealWeightInverseMeanRatio(err),
                               0 };

    // Read Mesh
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk("../../meshFiles/3D/VTK/12-pyramid-unit-sphere.vtk", err);
  CPPUNIT_ASSERT(!err);

    // Check that the mesh read correctly, and contains what is
    // expected later.

    // Expecting file to contain 12 pyramid elements constructed
    // from 15 vertices.
  size_t num_vtx, num_pyr, conn_len;
  mesh->get_all_sizes( num_vtx, num_pyr, conn_len, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(num_vtx == 15);
  CPPUNIT_ASSERT(num_pyr == 12);
  CPPUNIT_ASSERT(conn_len == 60);
  
    // Get mesh data
  Mesh::VertexHandle vert_array[15];
  Mesh::ElementHandle elem_array[12];
  size_t conn_offsets[13], conn_indices[60];
  EntityTopology type_array[12];
  mesh->get_all_mesh( vert_array, 15, elem_array, 12, 
                      conn_offsets, 13, conn_indices, 60, err );
  CPPUNIT_ASSERT(!err);
  mesh->elements_get_topologies( elem_array, type_array, 12, err );
  CPPUNIT_ASSERT(!err);
  
    // Verify element types and number of vertices
  for (unsigned i = 0; i < 12; ++i)
  {
    CPPUNIT_ASSERT( type_array[i] == PYRAMID );
    CPPUNIT_ASSERT( conn_offsets[i] == 5*i );
  }
  
    // All pyramids should share a common apex, at the
    // center of the sphere
  size_t apex_index = conn_indices[4];
  for (unsigned i = 1; i < 12; ++i)
  {
    CPPUNIT_ASSERT( conn_indices[5*i+4] == apex_index );
  }
  
    // Verify that apex is at origin and all other vertices are
    // on unit sphere
  MsqVertex vertices[15];
  mesh->vertices_get_coordinates( vert_array, vertices, 15, err );
  CPPUNIT_ASSERT(!err);
  for (unsigned i = 0; i < 15; ++i)
  {
    if (i == apex_index)
      CPPUNIT_ASSERT( vertices[i].within_tolerance_box( Vector3D(0,0,0), 1e-6 ) );
    else
      CPPUNIT_ASSERT( fabs(1.0 - vertices[i].length()) < 1e-6 );
  }
  
    // Try smoothing w/out moving the free vertex and verify that
    // the smoother didn't move the vertex
  Vector3D position(0,0,0);
  for (unsigned i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( mesh, vert_array[apex_index], position, metrics[i] ) );
  
    // Now try moving the vertex and see if the smoother moves it back
    // to the origin
  position.set( 0.5, 0.5, 0.5 );
  for (unsigned i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( mesh, vert_array[apex_index], position, metrics[i] ) );

  return 0;
}
  
  
bool smooth_mesh( Mesh* mesh,
                  Mesh::VertexHandle free_vertex_at_origin, 
                  Vector3D initial_free_vertex_position,
                  QualityMetric* metric )
{
  Mesquite::MsqPrintError err(cout);
  const Vector3D origin( 0, 0, 0 );
  
  // print a little output so we know when we died
  msq_stdio::cout << 
  "**************************************************************************" 
  << msq_stdio::endl << 
  "* Smoothing..."
  << msq_stdio::endl << 
  "* Metric: " << metric->get_name()
  << msq_stdio::endl << 
  "* Apex position: " << initial_free_vertex_position
  << msq_stdio::endl << 
  "**************************************************************************" 
  << msq_stdio::endl;
    
  // Use numeric approx of derivitives until analytic solutions
  // are working for pyramids
  //metric->set_gradient_type( QualityMetric::NUMERICAL_GRADIENT );
  //metric->set_hessian_type( QualityMetric::NUMERICAL_HESSIAN );
  
  
  // Set free vertex to specified position
  mesh->vertex_set_coordinates( free_vertex_at_origin, 
                                initial_free_vertex_position,
                                err );
  CPPUNIT_ASSERT(!err);

  // Create a MeshSet object and InstructionQueue
  MeshSet mesh_set;
  mesh_set.add_mesh(mesh, err); 
  CPPUNIT_ASSERT(!err);
  InstructionQueue Q;

  // Set up objective function
  LPtoPTemplate* obj_func = new LPtoPTemplate(metric, 1, err);
  CPPUNIT_ASSERT(!err);
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);

  // Create solver
  VertexMover* solver = new ConjugateGradient( obj_func, err );
  CPPUNIT_ASSERT(!err);
  solver->set_patch_type(PatchData::GLOBAL_PATCH, err);
  CPPUNIT_ASSERT(!err);

  // Set stoping criteria for solver
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_double( 
    TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, 1e-6, err);
  CPPUNIT_ASSERT(!err);
  solver->set_inner_termination_criterion(&tc_inner);
   
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  CPPUNIT_ASSERT(!err);
  solver->set_outer_termination_criterion(&tc_outer);

  // Need to do this too
  solver->add_culling_method(PatchData::NO_BOUNDARY_VTX);
   
  // Add solver to queue
  Q.set_master_quality_improver(solver, err); 
  CPPUNIT_ASSERT(!err);
 
  // And smooth...
  Q.run_instructions(mesh_set, err); 
  CPPUNIT_ASSERT(!err);
  
  // Verify that vertex was moved back to origin
  MsqVertex vtx;
  mesh->vertices_get_coordinates( &free_vertex_at_origin, &vtx, 1, err );
  CPPUNIT_ASSERT( !err );
  Vector3D position = vtx;
  
  // print a little output so we know when we died
  msq_stdio::cout << 
  "**************************************************************************" 
  << msq_stdio::endl << 
  "* Done Smoothing:"
  << msq_stdio::endl << 
  "* Metric: " << metric->get_name()
  << msq_stdio::endl << 
  "* Apex position: " << position
  << /*msq_stdio::endl << */ 
  "**************************************************************************" 
  << msq_stdio::endl;
  
  CPPUNIT_ASSERT( position.within_tolerance_box( Vector3D(0,0,0), 1e-6 ) );
  return false;
}
