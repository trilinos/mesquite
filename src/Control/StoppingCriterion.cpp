// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
// DESCRIPTION:
// ============
/*! \file StoppingCriterion.hpp
  
    \brief  Member function of the Mesquite::StoppingCriterion class

    \author Thomas Leurent
    \date   2002-07-23
 */



#include "StoppingCriterion.hpp"
#include "QualityAssessor.hpp"
using namespace Mesquite;
int StoppingCriterion::loopCounter;

#undef __FUNC__
#define  __FUNC__ "StoppingCriterion::StoppingCriterion"

StoppingCriterion::StoppingCriterion(ObjectiveFunction* of,
                                     double lower, double upper) :
   criterion(of),
   lowerBound(lower),
   upperBound(upper)
{

}


#undef __FUNC__
#define  __FUNC__ "StoppingCriterion::StoppingCriterion"

StoppingCriterion::StoppingCriterion(QualityAssessor* qm,
                                     double lower, double upper) :
   criterion(qm),
   lowerBound(lower),
   upperBound(upper)
{

}


#undef __FUNC__
#define  __FUNC__ "StoppingCriterion::StoppingCriterion"
StoppingCriterion::StoppingCriterion(enum SCFunction func,
                                     double lower, double upper):
    criterion(func),
   lowerBound(lower),
   upperBound(upper)
{
   
}

#undef __FUNC__
#define  __FUNC__ "StoppingCriterion::StoppingCriterion"
StoppingCriterion::StoppingCriterion(enum SCFunction func,
                                     int alpha):
    criterion(func),
    mAlpha(alpha)
{
  loopCounter=0;
}



bool StoppingCriterion::stop(MeshSet &ms, MsqError &err)
{
  double value=0;
   switch( criterion.mType ) {
   case StoppingCriterionEntry::OBJ_FUNC_CRIT:
        //double value = criterion.mObjFunc->evaluate_mesh( ms,depth,err );
        return true;
   case StoppingCriterionEntry::ASSESSOR_CRIT:
      value = criterion.mAssessor->assess_mesh_quality(ms,err);
      if(value<upperBound && value>lowerBound)
        return true;
      else
        return false;
      break;
   case StoppingCriterionEntry::SCFUNCTION_CRIT:
      switch(criterion.mSCFunc){
        case(NUMBER_OF_PASSES): 
           if(loopCounter>=mAlpha)
             return true;
           else
             return false;
           break;
        default:
           err.set_msg("Enum not implemented.");
           return true;
      };
      
        
        // return ...;
   default:
      err.set_msg("StoppingCriterion constructors must prevent from getting here.");
   }
   
   return false;
}

CompositeAndStoppingCriterion::CompositeAndStoppingCriterion(
   StoppingCriterion* SC1, StoppingCriterion* SC2 ) 
{
    stopCrit1=SC1;
    stopCrit2=SC2;
}

#undef __FUNC__
#define  __FUNC__ "CompositeAndStoppingCriterion::stop"
bool CompositeAndStoppingCriterion::stop(MeshSet &ms, MsqError &err)
{
   bool stop1 = stopCrit1->stop(ms, err); MSQ_CHKERR(err);
   if(!stop1)
     return false;
   bool stop2 = stopCrit2->stop(ms, err); MSQ_CHKERR(err);
   return (stop2);
}

CompositeOrStoppingCriterion::CompositeOrStoppingCriterion(
   StoppingCriterion* SC1, StoppingCriterion* SC2 ) :
         stopCrit1(SC1), stopCrit2(SC2) { }

#undef __FUNC__
#define  __FUNC__ "CompositeOrStoppingCriterion::stop"
bool CompositeOrStoppingCriterion::stop(MeshSet &ms, MsqError &err)
{
   bool stop1 = stopCrit1->stop(ms, err); MSQ_CHKERR(err);
   if(stop1)
     return stop1;
   stop1 = stopCrit2->stop(ms, err); MSQ_CHKERR(err);
   return (stop1);
}

