// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ShapeImprovementWrapper.cpp

Member functions of the Mesquite::ShapeImprovementWrapper class

  \author Michael Brewer
  \date   June 6, 2003
 */

#include "InstructionQueue.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "ShapeImprovementWrapper::ShapeImprovementWrapper"
/*! The consturctor allows for two values.  The first is a 
  time bound (in seconds) used as a termination criterion.  If
  this value is non-positive, no time bound will be set.
  By default, the value is set to zero and no time bound
  is used.  The second value is the tolerance for the gradient
  norm termination criteria.  The default value is 1.e-6.*/
ShapeImprovementWrapper::ShapeImprovementWrapper(double cpu_time,
                                                 double grad_norm) {

    //arbitraryily chosen variables
  untBeta=1.e-8;
  successiveEps=1.e-1;
  
  
  
  if(cpu_time>0.0){
    timerNeeded=true;
  }
  else{
    timerNeeded=false;
  }
  maxTime=cpu_time;
  
  MsqError err;
  untangleMetric = new UntangleBetaQualityMetric(untBeta);
  untangleFunc =  new LPtoPTemplate(untangleMetric, 2, err);
  untangleFunc->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  untangleGlobal = new ConjugateGradient(untangleFunc,err);
  untangleGlobal->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
  
  untangleLocal = new ConjugateGradient(untangleFunc,err);
  untangleLocal->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err,1 ,1);
  
  untangleGlobalInner = new TerminationCriterion();
  untangleGlobalOuter = new TerminationCriterion();
  untangleLocalInner = new TerminationCriterion();
  untangleLocalOuter = new TerminationCriterion();
  
  untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE,0.0,err);
  untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE,successiveEps,err);
  untangleGlobalOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  
  untangleLocalInner->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,3,err);
  untangleLocalInner->set_culling_type(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, 0.0 , err);
  
  untangleLocalOuter->add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE,0.0,err);
  untangleLocalOuter->add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE,successiveEps*.1,err);
  
  
  meanRatio = new MeanRatioQualityMetric;
  meanRatio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
  meanRatio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
  meanRatio->set_averaging_method(QualityMetric::LINEAR,err);
    // creates the l_2 squared objective function
  objFunc = new LPtoPTemplate(meanRatio, 2, err);
  objFunc->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
    //creates a FeasibleNewtone improver
  feasNewt = new FeasibleNewton(objFunc);
  feasNewt->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
  mQA = new QualityAssessor(meanRatio,QualityAssessor::MAXIMUM);
  mQA->add_quality_assessment(meanRatio, QualityAssessor::MINIMUM,err);
  mQA->add_quality_assessment(meanRatio, QualityAssessor::AVERAGE,err);
  mQA->add_quality_assessment(meanRatio, QualityAssessor::RMS,err);   
        //**************Set stopping criterion*e***************
  termInner = new TerminationCriterion();
  termOuter = new TerminationCriterion();
  termInner->add_criterion_type_with_double(TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE,grad_norm,err);
  termOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
    // sets a culling method on the first QualityImprover
  untangleGlobal->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  untangleGlobal->set_inner_termination_criterion(untangleGlobalInner);
  untangleGlobal->set_outer_termination_criterion(untangleGlobalOuter);
    // sets a culling method on the second QualityImprover
  untangleLocal->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  untangleLocal->set_inner_termination_criterion(untangleLocalInner);
  untangleLocal->set_outer_termination_criterion(untangleLocalOuter);
    // sets a culling method on the third QualityImprover
  feasNewt->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  feasNewt->set_inner_termination_criterion(termInner);
  feasNewt->set_outer_termination_criterion(termOuter);
      
}


#undef __FUNC__
#define __FUNC__ "ShapeImprovementWrapper::~ShapeImprovementWrapper"
ShapeImprovementWrapper::~ShapeImprovementWrapper()
{
  delete untangleMetric;
  delete untangleFunc;
  delete untangleGlobal;
  delete untangleLocal;
  delete untangleGlobalInner;
  delete untangleGlobalOuter;
  delete untangleLocalInner;
  delete untangleLocalOuter;
      
  delete meanRatio;
  delete objFunc;
  delete feasNewt;
  delete mQA;
  delete termInner;
  delete termOuter;
}


#undef __FUNC__
#define __FUNC__ "ShapeImprovementWrapper::run_instructions"
/*!Run instructions first calls the global untangler.  If the
  resulting mesh is tangled after that pre-conditioning step,
  The mesh is iteratively smoothed with a local and then global
  untangler until the mesh is untangled or until a certain time
  constraint has been exceeded.  If the mesh was successfully
  untangled and there is still time remaining, a mean ratio
  shape improvement is then performed.*/
void ShapeImprovementWrapper::run_instructions(MeshSet &ms, MsqError &err)
{
    //a timer to keep track of the amount of time spent in this wrapper
  Timer totalTimer;
    //time remaining keeps a track of how much time is left before the
    //wrapper must terminate.  If the wrapper is set to terminate on
    //a time constraint, time_remaining will always be 1.0
  double time_remaining=1.0;
  mQA->loop_over_mesh(ms, err);
    //if using a time constraint set the termination criteria.
  if(timerNeeded){
    time_remaining=maxTime;
    untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
  }
    //global untangler
  untangleGlobal->loop_over_mesh(ms, err);
  if(timerNeeded)
    time_remaining=maxTime-totalTimer.since_birth();
  double func_val=untangleGlobalInner->get_current_function_value();
    //if there is time remaining and the mesh is tangled iterate
    //over untanglers
  while(func_val> 0.0 && time_remaining>0){
    if(timerNeeded)
      untangleLocalOuter->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
    untangleLocal->loop_over_mesh(ms, err);
    func_val=untangleGlobalInner->get_current_function_value();
    if(timerNeeded)
      time_remaining=maxTime-totalTimer.since_birth();
    if(func_val>0 && time_remaining>0){
      if(timerNeeded)
        untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
      untangleGlobal->loop_over_mesh(ms, err);
      func_val=untangleGlobalInner->get_current_function_value();
      if(timerNeeded)
        time_remaining=maxTime-totalTimer.since_birth();
    }
  }
  mQA->loop_over_mesh(ms, err);
  if(timerNeeded)
    time_remaining=maxTime-totalTimer.since_birth();
    //if all the time constraint has been exceeded, notify the user that
    //the shape improvement has not been performed.
  if(time_remaining<=0){
    PRINT_INFO("\nOptimization is terminating without perfoming shape improvement");
    PRINT_INFO("\n Untangle Function Value is %f",func_val);
  }
    //otherwise, perform the shape improvement.
  else{
    if(timerNeeded)
      termInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
    feasNewt->loop_over_mesh(ms, err);
    mQA->loop_over_mesh(ms, err);
  } 
}

  
