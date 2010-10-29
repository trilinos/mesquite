/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ShapeImprover.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */


#include "ShapeImprover.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "PMeanPTemplate.hpp"
#include "ConjugateGradient.hpp"
#include "TerminationCriterion.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"

#include "TRel2DShapeBarrier.hpp"
#include "TRel2DShape.hpp"
#include "TRel3DShapeBarrier.hpp"
#include "TRel3DShape.hpp"
#include "TRelQualityMetric.hpp"
#include "IdealShapeTarget.hpp"

namespace MESQUITE_NS {

const double DEFAULT_BETA = 0.005;
const int DEFUALT_PARALLEL_ITERATIONS = 10;

ShapeImprover::ShapeImprover() 
 : maxTime(-1.0), 
   mBeta(DEFAULT_BETA),
   parallelIterations(DEFUALT_PARALLEL_ITERATIONS)
{}

void ShapeImprover::set_vertex_movement_limit_factor( double beta )
{
  assert(beta > 0.0);
  assert(beta < 1.0);
  mBeta = beta;
}

void ShapeImprover::set_cpu_time_limit( double limit )
{
  assert(limit >= 0.0);
  maxTime = limit;
}

void ShapeImprover::set_parallel_iterations( int count )
{
  if (count < 1) {
    assert(false);
    count = 1;
  }
  parallelIterations = count;
}

void ShapeImprover::run_wrapper( Mesh* mesh,
                                ParallelMesh* pmesh,
                                MeshDomain* domain,
                                Settings* settings,
                                QualityAssessor* qa,
                                MsqError& err )
{
    // Quality Metrics
  IdealShapeTarget target;
    // No-barrier phase
  TRel2DShape mu_2_no;
  TRel3DShape mu_3_no;
  TRelQualityMetric metric_no( &target, &mu_2_no, &mu_3_no );
    // Barrier phase
  TRel2DShapeBarrier mu_2_b;
  TRel3DShapeBarrier mu_3_b;
  TRelQualityMetric metric_b( &target, &mu_2_b, &mu_3_b );
    // QualityAssessor
  qa->add_quality_assessment( &metric_no );
  qa->add_quality_assessment( &metric_b );
  

    // Phase 1: use non-barrier metric until mesh does not
    //          have any inverted elements
  PMeanPTemplate obj_func_no( 1.0, &metric_no );
  ConjugateGradient improver_no( &obj_func_no );
  improver_no.use_global_patch();
  TerminationCriterion inner_no;
  inner_no.add_absolute_vertex_movement_edge_length( mBeta );
  inner_no.add_untangled_mesh();
  if (maxTime > 0.0)
    inner_no.add_cpu_time( maxTime );
  improver_no.set_inner_termination_criterion( &inner_no );
  InstructionQueue q_no;
  Timer totalTimer;
  q_no.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q_no.set_master_quality_improver( &improver_no, err ); MSQ_ERRRTN(err);
  q_no.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q_no.run_common( mesh, pmesh, domain, settings, err ); MSQ_ERRRTN(err);
  
    // Phase 2: use-barrer metric on untangled mesh
  PMeanPTemplate obj_func_b( 1.0, &metric_b );
  ConjugateGradient improver_b( &obj_func_b );
  improver_b.use_global_patch();
  TerminationCriterion inner_b;
  inner_b.add_absolute_vertex_movement_edge_length( mBeta );
  if (maxTime > 0.0) {
    double remaining = maxTime - totalTimer.since_birth();
    if (remaining <= 0.0 ){
      //MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement." << std::endl;
      remaining = 0.0;
    }
    inner_b.add_cpu_time( remaining );
  }
  improver_b.set_inner_termination_criterion( &inner_b );
  InstructionQueue q_b;
  q_b.set_master_quality_improver( &improver_b, err ); MSQ_ERRRTN(err);
  q_b.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q_b.run_common( mesh, pmesh, domain, settings, err ); MSQ_ERRRTN(err);
}

} // namespace MESQUITE_NS
