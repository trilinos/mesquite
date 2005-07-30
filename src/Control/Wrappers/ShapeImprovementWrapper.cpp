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
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ShapeImprovementWrapper.cpp

Member functions of the Mesquite::ShapeImprovementWrapper class

  \author Michael Brewer
  \date   June 6, 2003
 */

#include "InstructionQueue.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

namespace Mesquite {

/*! The consturctor allows for two values.  The first is a 
  time bound (in seconds) used as a termination criterion.  If
  this value is non-positive, no time bound will be set.
  By default, the value is set to zero and no time bound
  is used.  The second value is the tolerance for the gradient
  norm termination criteria.  The default value is 1.e-6.*/
ShapeImprovementWrapper::ShapeImprovementWrapper(MsqError& err,
                                                 double cpu_time,
                                                 double grad_norm) 
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
  
  
  
  if(cpu_time>0.0){
    timerNeeded=true;
  }
  else{
    timerNeeded=false;
  }
  maxTime=cpu_time;
  
  untangleMetric = new UntangleBetaQualityMetric(untBeta);
  untangleFunc =  new LPtoPTemplate(untangleMetric, 2, err);  MSQ_ERRRTN(err);
  untangleFunc->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  untangleGlobal = new ConjugateGradient(untangleFunc,err);  MSQ_ERRRTN(err);
  untangleGlobal->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);  MSQ_ERRRTN(err);
  
  untangleGlobalInner = new TerminationCriterion();
  untangleGlobalOuter = new TerminationCriterion();
  
  untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE,0.0,err);  MSQ_ERRRTN(err);
  untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_ABSOLUTE,successiveEps,err);  MSQ_ERRRTN(err);
  untangleGlobalOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);  MSQ_ERRRTN(err);
  
  inverseMeanRatio = new IdealWeightInverseMeanRatio(err); MSQ_ERRRTN(err);
  inverseMeanRatio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
  inverseMeanRatio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
  inverseMeanRatio->set_averaging_method(QualityMetric::LINEAR,err);  MSQ_ERRRTN(err);
    // creates the l_2 squared objective function
  objFunc = new LPtoPTemplate(inverseMeanRatio, 2, err);  MSQ_ERRRTN(err);
  objFunc->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
    //creates a FeasibleNewtone improver
  feasNewt = new FeasibleNewton(objFunc);
  feasNewt->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);  MSQ_ERRRTN(err);
  mQA = new QualityAssessor(inverseMeanRatio,QualityAssessor::MAXIMUM, err); MSQ_ERRRTN(err);
  mQA->add_quality_assessment(inverseMeanRatio, QualityAssessor::MINIMUM,err);  MSQ_ERRRTN(err);
  mQA->add_quality_assessment(inverseMeanRatio, QualityAssessor::AVERAGE,err);  MSQ_ERRRTN(err);
  mQA->add_quality_assessment(inverseMeanRatio, QualityAssessor::RMS,err);  MSQ_ERRRTN(err);   
        //**************Set stopping criterion*e***************
  termInner = new TerminationCriterion();
  termOuter = new TerminationCriterion();
  termInner->add_criterion_type_with_double(TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE,grad_norm,err);  MSQ_ERRRTN(err);
  termInner->add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE,successiveEps,err);  MSQ_ERRRTN(err);
  termOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);  MSQ_ERRRTN(err);
  untangleGlobal->set_inner_termination_criterion(untangleGlobalInner);
  untangleGlobal->set_outer_termination_criterion(untangleGlobalOuter);
    //untangleLocal->set_inner_termination_criterion(untangleLocalInner);
    //untangleLocal->set_outer_termination_criterion(untangleLocalOuter);
  feasNewt->set_inner_termination_criterion(termInner);
  feasNewt->set_outer_termination_criterion(termOuter);
      
}


ShapeImprovementWrapper::~ShapeImprovementWrapper()
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
  untangler until the mesh is untangled or until a certain time
  constraint has been exceeded.  If the mesh was successfully
  untangled and there is still time remaining, an inverse mean ratio
  shape improvement is then performed.*/
void ShapeImprovementWrapper::run_instructions( Mesh* mesh,
                                                MeshDomain* domain, 
                                                MsqError &err)
{
    //a timer to keep track of the amount of time spent in this wrapper
  Timer totalTimer;
    //time remaining keeps a track of how much time is left before the
    //wrapper must terminate.  If the wrapper is set to terminate on
    //a time constraint, time_remaining will always be 1.0
  double time_remaining=1.0;
  mQA->loop_over_mesh(mesh, domain, 0, err);  MSQ_ERRRTN(err);
    //if using a time constraint set the termination criteria.
  if(timerNeeded){
    time_remaining=maxTime;
    untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
    MSQ_ERRRTN(err);
  }
    //global untangler
  untangleGlobal->loop_over_mesh(mesh, domain, 0, err);  MSQ_ERRRTN(err);
  if(timerNeeded)
    time_remaining=maxTime-totalTimer.since_birth();
  double func_val=untangleGlobalInner->get_current_function_value();

  mQA->loop_over_mesh(mesh, domain, 0, err);  MSQ_ERRRTN(err);
  if(timerNeeded)
    time_remaining=maxTime-totalTimer.since_birth();
    //if all the time constraint has been exceeded, notify the user that
    //the shape improvement has not been performed.
  if(time_remaining<=0){
    MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement.  "
                  << "Untangle Function Value is " << func_val << " .\n";
  }
    //otherwise, perform the shape improvement.
  else{
    if(timerNeeded) {
      termInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
      MSQ_ERRRTN(err);
    }
    feasNewt->loop_over_mesh( mesh, domain, 0, err);  MSQ_ERRRTN(err);
    mQA->loop_over_mesh(mesh, domain, 0, err);  MSQ_ERRRTN(err);
  } 
}

} // namespace Mesquite

  
