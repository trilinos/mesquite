/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ViscousCFDTetShapeWrapper.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ViscousCFDTetShapeWrapper.hpp"

#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TagVertexMesh.hpp"
#include "TrustRegion.hpp"
#include "TerminationCriterion.hpp"

#include "PMeanPTemplate.hpp"
#include "TMPQualityMetric.hpp"
#include "AddQualityMetric.hpp"

#include "Target3DShapeBarrier.hpp"
#include "Target2DShapeBarrier.hpp"
#include "Target3DShape.hpp"
#include "Target2DShape.hpp"
#include "Target3DShapeSizeOrientBarrier.hpp"
#include "Target2DShapeSizeOrientBarrier.hpp"
#include "Target3DShapeSizeOrient.hpp"
#include "Target2DShapeSizeOrient.hpp"

#include "IdealTargetCalculator.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "ReferenceMesh.hpp"
#include "TetDihedralWeight.hpp"
#include "RemainingWeight.hpp"

namespace MESQUITE_NS {

void ViscousCFDTetShapeWrapper::run_instructions( Mesh* mesh, 
                                                  MeshDomain* domain, 
                                                  MsqError& err )
{
  InstructionQueue q( *this ); // copy settings to queue
  
  // Set up barrier metric to see if mesh contains inverted elements
  Target3DShapeBarrier mu3Db;
  Target2DShapeBarrier mu2Db;
  IdealTargetCalculator w_ideal;
  TMPQualityMetric barrier( &w_ideal, &mu2Db, &mu3Db );
  
  // Check for inverted elements in the mesh
  QualityAssessor inv_check( &barrier );
  inv_check.disable_printing_results();
  q.add_quality_assessor( &inv_check, err );  MSQ_ERRRTN(err);
  q.run_instructions( mesh, domain, err ); MSQ_ERRRTN(err);
  q.remove_quality_assessor( 0, err ); MSQ_ERRRTN(err);
  const QualityAssessor::Assessor* inv_b = inv_check.get_results( &barrier );
  const bool use_barrier = (0 == inv_b->get_invalid_element_count());
  
  // Create remaining metric instances
  Target3DShape mu3D;
  Target2DShape mu2D;
  Target3DShapeSizeOrient mu3Do;
  Target2DShapeSizeOrient mu2Do;
  Target3DShapeSizeOrientBarrier mu3Dob;
  Target2DShapeSizeOrientBarrier mu2Dob;
  
  // Select which target metrics to use
  TargetMetric3D *mu3Dp, *mu3Dop;
  TargetMetric2D *mu2Dp, *mu2Dop;
  if (use_barrier) {
    mu3Dp = &mu3Db;
    mu2Dp = &mu2Db;
    mu3Dop = &mu3Dob;
    mu2Dop = &mu2Dob;
  }
  else {
    mu3Dp = &mu3D;
    mu2Dp = &mu2D;
    mu3Dop = &mu3Do;
    mu2Dop = &mu2Do;
  }
  
  // Set up target and weight calculators
  TagVertexMesh init_mesh( err, mesh );  MSQ_ERRRTN(err);
  ReferenceMesh ref_mesh( &init_mesh );
  RefMeshTargetCalculator w_init( &ref_mesh );
  TetDihedralWeight c_dihedral( &ref_mesh, dCutoff, aVal );
  RemainingWeight c_remaining( &c_dihedral );
  
  // Create objective function
  TMPQualityMetric metric1( &w_ideal, &c_dihedral,  mu2Dp,  mu3Dp  );
  TMPQualityMetric metric2( &w_init,  &c_remaining, mu2Dop, mu3Dop );
  AddQualityMetric of_metric( &metric2, &metric2, err );  MSQ_ERRRTN(err);
  PMeanPTemplate obj_func( 1.0, &of_metric );
  
  // Create optimizer
  TrustRegion solver( &obj_func );
  TerminationCriterion term;
  term.add_iteration_limit( iterationLimit );
  term.add_absolute_vertex_movement( maxVtxMovement );
  solver.set_inner_termination_criterion( &term );
  
  // Create instruction queue
  QualityAssessor qa;
  qa.add_quality_assessment( &metric1 );
  qa.add_quality_assessment( &metric2 );
  qa.add_quality_assessment( &of_metric );
  q.add_quality_assessor( &qa, err ); MSQ_ERRRTN(err);
  q.set_master_quality_improver( &solver, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( &qa, err ); MSQ_ERRRTN(err);

  // Optimize mesh
  q.run_instructions( mesh, domain, err ); MSQ_CHKERR(err);  
}

} // namespace MESQUITE_NS
