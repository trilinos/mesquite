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
//  LAST-MOD: 14-Nov-02 at 17:31:53 by Thomas Leurent


/*! \file LaplacianIQ.cpp

This is the second possibility for wrappers. It is based on the InctructionQueue concept. 

 */
// DESCRIP-END.
//


#ifndef LaplacianIQ_hpp
#define LaplacianIQ_hpp

#include "MeanRatioQualityMetric.hpp" 
#include "LaplacianSmoother.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "StoppingCriterion.hpp"

namespace Mesquite { 

   class LaplacianIQ : public InstructionQueue {
   private:
      ShapeQualityMetric* meanRatio;
      LaplacianSmoother* lapl1;
      QualityAssessor* mQA;
      StoppingCriterion* mStop;

   public:
      
      //! Constructor sets the instructions in the queue.  
      LaplacianIQ() {
         MsqError err;
         // creates a mean ratio quality metric ...
         meanRatio = MeanRatioQualityMetric::create_new();
     
         // creates the laplacian smoother  procedures
         lapl1 = new LaplacianSmoother();
         mQA = new QualityAssessor(meanRatio,QualityAssessor::MAXIMUM);
     
         //**************Set stopping criterion****************
            mStop = new StoppingCriterion(StoppingCriterion::NUMBER_OF_PASSES,10);
            lapl1->set_stopping_criterion(mStop);
            // sets a culling method on the first QualityImprover
            lapl1->add_culling_method(QualityImprover::NO_BOUNDARY_VTX);
      
            // adds 1 pass of pass1 
            this->add_quality_assessor(mQA,err); MSQ_CHKERR(err);
            this->set_master_quality_improver(lapl1, err); MSQ_CHKERR(err);
            this->add_quality_assessor(mQA,err); MSQ_CHKERR(err);
      }

      
      //! Destructor must delete the objects inserted in the queue.
      ~LaplacianIQ()
      {
         delete meanRatio;
         delete lapl1;
         delete mQA;
         delete mStop;
      }
  
   };


} // namespace

#endif // LaplacianIQ_hpp
