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

#include "UntangleBetaQualityMetric.hpp"
#include "ConjugateGradient.hpp"

namespace Mesquite { 
  /*! \class ShapeImprovementWrapper
       \brief Wrapper which performs a Feasible Newton solve using
       an \f$\ell_2^2 \f$ objective function template with mean
       ratio.
       
     */
  class ShapeImprovementWrapper : public InstructionQueue {
     
  public:  
      //Constructor sets the instructions in the queue.
    ShapeImprovementWrapper(double cpu_time = 0.0, double grad_norm =1.e-6);
    
      //! Destructor must delete the objects inserted in the queue.
    virtual ~ShapeImprovementWrapper();
    
      //! run_instructions runs the wrapper on the given MeshSet.
    virtual void run_instructions(MeshSet &ms, MsqError &err);
    
    
  private:
    UntangleQualityMetric* untangleMetric;
    LPtoPTemplate* untangleFunc;
    
    VertexMover* untangleGlobal;
    
    TerminationCriterion* untangleGlobalOuter; 
    TerminationCriterion* untangleGlobalInner;

    ShapeQualityMetric* meanRatio; 
    LPtoPTemplate* objFunc;
    FeasibleNewton* feasNewt;
    QualityAssessor* mQA;
    TerminationCriterion* termOuter; 
    TerminationCriterion* termInner;

    bool timerNeeded;
    double maxTime;
      //arbitraryily chosen variables
    double untBeta;
    double successiveEps;
  

  };
  
  
} // namespace

#endif // ShapeImprovementWrapper_hpp
