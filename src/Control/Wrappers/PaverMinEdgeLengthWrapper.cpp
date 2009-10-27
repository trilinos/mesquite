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


/** \file PaverMinEdgeLengthWrapper.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "PaverMinEdgeLengthWrapper.hpp"

#include "TerminationCriterion.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"
#include "MeshImpl.hpp"
#include "PlanarDomain.hpp"

#include "PMeanPTemplate.hpp"
#include "TrustRegion.hpp"
#include "TMPQualityMetric.hpp"
#include "IdealTargetCalculator.hpp"
#include "Target2DShapeSizeBarrier.hpp"
//#include "Target3DShapeSizeBarrier.hpp"
#include "TargetSize.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "ReferenceMesh.hpp"

#include "EdgeLengthMetric.hpp"
#include "LambdaConstant.hpp"

namespace MESQUITE_NS {

static double calculate_average_lambda( Mesh* mesh, 
                                        TargetCalculator* tc,
                                        MsqError& err );

void PaverMinEdgeLengthWrapper::run_instructions_internal( Mesh* mesh, 
                                                  ParallelMesh* pmesh,
                                                  MeshDomain* domain, 
                                                  MsqError& err )
{
  InstructionQueue q( *this ); // copy settings to queue
 
    // calculate average lambda for mesh
  ReferenceMesh ref_mesh( mesh );
  RefMeshTargetCalculator W_0( &ref_mesh );
  double lambda = calculate_average_lambda( mesh, &W_0, err ); MSQ_ERRRTN(err);
  
    // create objective function
  IdealTargetCalculator W_i;
  LambdaConstant W( lambda, &W_i );
  Target2DShapeSizeBarrier tm2;
  //Target3DShapeSizeBarrier tm3;
  TMPQualityMetric mu( &W, &tm2, 0 /*&tm3*/ );
  PMeanPTemplate of( 1.0, &mu );
  
    // create quality assessor
  EdgeLengthMetric len(0.0);
  QualityAssessor qa( &mu );
  qa.add_quality_assessment( &len );
  q.add_quality_assessor( &qa, err );
  
    // create solver
  TrustRegion solver( &of );
  TerminationCriterion tc;
  tc.add_absolute_vertex_movement( 1e-2 );
  solver.set_inner_termination_criterion( &tc );
  q.set_master_quality_improver( &solver, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( &qa, err );

  // Optimize mesh
  if (pmesh)
    q.run_instructions( pmesh, domain, err ); 
  else
    q.run_instructions( mesh, domain, err ); 
  MSQ_CHKERR(err);  
}

double calculate_average_lambda( Mesh* mesh, 
                                 TargetCalculator* tc,
                                 MsqError& err )
{
  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err );
  
  std::vector<size_t> handles;
  TMPQualityMetric::get_patch_evaluations( pd, handles, false, err );

  double lambda = 0.0;
  for (size_t i = 0; i < handles.size(); ++i) {
    const Sample s = ElemSampleQM::sample( handles[i] );
    const size_t e = ElemSampleQM::  elem( handles[i] );
    if (TopologyInfo::dimension( pd.element_by_index(e).get_element_type() ) == 2)
    {
      MsqMatrix<3,2> W;
      tc->get_2D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
      lambda += TargetSize::factor_size( W );
    }
    else
    {
      MsqMatrix<3,3> W;
      tc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
      lambda += TargetSize::factor_size( W );
    }
  }
  return lambda/handles.size();
}

} // namespace MESQUITE_NS
