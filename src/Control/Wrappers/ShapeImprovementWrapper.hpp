// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 14-Nov-02 at 16:51:36
//  LAST-MOD: 22-May-03 at 09:04:07 by Michael Brewer


/*! \file ShapeImprovementWrapper.hpp
  ShapeImprovementWrapper header file.

*/
// DESCRIP-END.
//


#ifndef ShapeImprovementWrapper_hpp
#define ShapeImprovementWrapper_hpp

#include "MeanRatioQualityMetric.hpp" 
#include "FeasibleNewton.hpp"
#include "LPtoPTemplate.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"

namespace Mesquite { 
  /*! \class ShapeImprovementWrapper
       \brief Wrapper which performs a Feasible Newton solve using
       an \f$\ell_2^2 \f$ objective function template with mean
       ratio.
       

       \todo MB:  We should probably add a conditional preconditioner
       such that if the initial mesh is tangled we attempt to untangle
       it.
     */
  class ShapeImprovementWrapper : public InstructionQueue {
     
  private:
    ShapeQualityMetric* meanRatio;
     
    LPtoPTemplate* objFunc;
    FeasibleNewton* feasNewt;
    QualityAssessor* mQA;
    TerminationCriterion* termOuter; 
    TerminationCriterion* termInner;

  public:
      
    /*! \brief Constructor sets the instructions in the queue.

      The consturctor allows for two values.  The first is a 
      time bound (in seconds) used as a termination criterion.  If
      this value is non-positive, no time bound will be set.
      By default, the value is set to zero and no time bound
      is used.  The second value is the tolerance for the gradient
      norm termination criteria.  The default value is 1.e-6.*/
    ShapeImprovementWrapper(double cpu_time = 0.0, double grad_norm =1.e-6) {
      MsqError err;
      meanRatio = MeanRatioQualityMetric::create_new();
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
      mQA->add_quality_assessment(meanRatio, QualityAssessor::AVERAGE,err);
      
      //**************Set stopping criterion****************
	termInner = new TerminationCriterion();
	termOuter = new TerminationCriterion();
	if(cpu_time> 0.0)
	  termInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,cpu_time,err);
	termInner->add_criterion_type_with_double(TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE,grad_norm,err);
	termOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
	// sets a culling method on the first QualityImprover
	feasNewt->add_culling_method(PatchData::NO_BOUNDARY_VTX);
	feasNewt->set_inner_termination_criterion(termInner);
	feasNewt->set_outer_termination_criterion(termOuter);
	 
	// adds 1 pass of pass1 
	this->add_quality_assessor(mQA,err); MSQ_CHKERR(err);
	this->set_master_quality_improver(feasNewt, err); MSQ_CHKERR(err);
	this->add_quality_assessor(mQA,err); MSQ_CHKERR(err);
    }
    
    
    //! Destructor must delete the objects inserted in the queue.
    virtual ~ShapeImprovementWrapper()
    {
      delete meanRatio;
      delete objFunc;
      delete feasNewt;
      delete mQA;
      delete termInner;
      delete termOuter;
    }
    
  };
  
  
} // namespace

#endif // ShapeImprovementWrapper_hpp
