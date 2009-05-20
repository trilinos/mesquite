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
 
    isenburg@llnl.gov, diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ParallelShapeImprovementWrapper.cpp

Member functions of the Mesquite::ShapeImprovementWrapper class

  \author Martin Isenburg
  \date   May 14, 2009
 */

#include "ParallelShapeImprovementWrapper.hpp"

namespace MESQUITE_NS {

/*! The value grad_norm is the tolerance for the gradient
  norm termination criteria.  The default value is 1.e-6.*/
ParallelShapeImprovementWrapper::ParallelShapeImprovementWrapper(MsqError& err, double grad_norm) 
 : untangleMetric(0),
   untangleFunc(0),
   untangleGlobal(0),
   untangleGlobalOuter(0),
   untangleGlobalInner(0),
   inverseMeanRatio(0),
   objFunc(0),
   feasNewt(0),
   mQA(0),
   termOuter(0),
   termInner(0)
{
    //arbitrarily chosen variables
  untBeta=1.e-8;
  successiveEps=1.e-4;
  
  untangleMetric = new UntangleBetaQualityMetric(untBeta);
  untangleFunc =  new LPtoPTemplate(untangleMetric, 2, err);  MSQ_ERRRTN(err);
  untangleGlobal = new ConjugateGradient(untangleFunc,err);  MSQ_ERRRTN(err);
  untangleGlobal->use_global_patch();
  untangleGlobalInner = new TerminationCriterion();
  untangleGlobalOuter = new TerminationCriterion();  
  untangleGlobalInner->add_absolute_quality_improvement( 0.0 );
  untangleGlobalInner->add_absolute_successive_improvement( successiveEps );
  untangleGlobalOuter->add_iteration_limit( 1 );
  untangleGlobal->set_inner_termination_criterion(untangleGlobalInner);
  untangleGlobal->set_outer_termination_criterion(untangleGlobalOuter);
  
  inverseMeanRatio = new IdealWeightInverseMeanRatio(err); MSQ_ERRRTN(err);
  inverseMeanRatio->set_averaging_method(QualityMetric::LINEAR);
  mQA = new QualityAssessor(inverseMeanRatio);
  objFunc = new LPtoPTemplate(inverseMeanRatio, 2, err);  MSQ_ERRRTN(err);
  feasNewt = new FeasibleNewton(objFunc,true);
  feasNewt->use_global_patch();
  // feasNewt->use_element_on_vertex_patch();
  termInner = new TerminationCriterion();
  termOuter = new TerminationCriterion();
  // termInner->add_iteration_limit( 5 );
  termInner->add_absolute_gradient_L2_norm( grad_norm );
  termInner->add_relative_successive_improvement( successiveEps );
  termOuter->add_iteration_limit( 1 );
  feasNewt->set_inner_termination_criterion(termInner);
  feasNewt->set_outer_termination_criterion(termOuter);      
}

ParallelShapeImprovementWrapper::~ParallelShapeImprovementWrapper()
{
  delete untangleMetric;
  delete untangleFunc;
  delete untangleGlobal;
  delete untangleGlobalInner;
  delete untangleGlobalOuter;
  delete inverseMeanRatio;
  delete objFunc;
  delete feasNewt;
  delete mQA;
  delete termInner;
  delete termOuter;
}

/*!Run instructions first calls the global untangler.  If the
  resulting mesh is tangled after that pre-conditioning step,
  The mesh is iteratively smoothed with a local and then global
  untangler until the mesh is untangled.  If the mesh was successfully
  untangled an inverse mean ratio shape improvement is then performed.*/
void ParallelShapeImprovementWrapper::run_instructions( ParallelMesh* parallel_mesh,
							MeshDomain* domain, 
							MsqError &err)
{
  mQA->loop_over_mesh(parallel_mesh, domain, 0, err);  MSQ_ERRRTN(err);  
  untangleGlobal->loop_over_mesh(parallel_mesh, domain, 0, err);  MSQ_ERRRTN(err);
  mQA->loop_over_mesh(parallel_mesh, domain, 0, err);  MSQ_ERRRTN(err);
  feasNewt->loop_over_mesh( parallel_mesh, domain, 0, err);  MSQ_ERRRTN(err);
  mQA->loop_over_mesh(parallel_mesh, domain, 0, err);  MSQ_ERRRTN(err);
}

} // namespace Mesquite
