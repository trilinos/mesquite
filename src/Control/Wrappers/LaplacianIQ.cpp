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


/** \file LaplacianIQ.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LaplacianIQ.hpp"
#include "IdealWeightInverseMeanRatio.hpp" 
#include "LaplacianSmoother.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"

void MESQUITE_NS::LaplacianIQ::run_wrapper( Mesh* mesh,
                                            ParallelMesh* pmesh,
                                            MeshDomain* domain,
                                            Settings* settings,
                                            QualityAssessor* qa,
                                            MsqError& err )
{
  if (iterationLimit < 1) {
    MSQ_SETERR(err)(MsqError::INVALID_ARG, "Invalid iteration limit: %d", iterationLimit );
    return;
  }

  IdealWeightInverseMeanRatio qa_metric;
  qa->add_quality_assessment( &qa_metric );
  
  LaplacianSmoother smoother;
  TerminationCriterion outer_term;
  outer_term.add_iteration_limit( iterationLimit );
  smoother.set_outer_termination_criterion( &outer_term );
  
  InstructionQueue q;
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.set_master_quality_improver( &smoother, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.run_common( mesh, pmesh, domain, settings, err ); MSQ_ERRRTN(err);
}
