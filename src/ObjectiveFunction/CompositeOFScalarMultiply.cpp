/*!
  \file    CompositeOFScalarMultiply.cpp
  \brief  

  This Objective Function combines two Objective Functions by mulitplication
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "ObjectiveFunction.hpp"
#include "CompositeOFScalarMultiply.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::CompositeOFScalarMultiply"
/*!
Sets the QualityMetric pointer to the metric associated with Obj.  If
Obj requires a feasible region, then so does the new CompositeOFScalarMultiply
ObjectiveFunction.  However, if alp is less than zero, the new
ObjectiveFunction's negateFlag is the opposite of Obj's.  
  \param alp (double)
  \param Obj (ObjectiveFunction*)
 */
CompositeOFScalarMultiply::CompositeOFScalarMultiply(double alp, ObjectiveFunction* Obj){

  set_quality_metric(Obj->get_quality_metric());
  set_feasible(Obj->get_feasible_constraint());
  objFunc=Obj;
  alpha=alp;
  if(alp<0)
    set_negate_flag(-1);
  else if (alp>0)
    set_negate_flag(1);
  else
    PRINT_WARNING("ObjectiveFunction being scaled by zero.");
}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::~CompositeOFScalarMultiply"

//Michael:  need to clean up here
CompositeOFScalarMultiply::~CompositeOFScalarMultiply(){

}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::concrete_evaluate"
/*!Computes fval= alpha*objFunc->evaluate(patch,err).  Note that since Obj's
  evaluate() function is called (as opposed to its concrete_evaluate) the
  returned value has been multiplied by objFunc's negateFlag (that is,
  if objFunc needed to be maximized then the value has been multiplied
  by negative one so that it may be minimized instead.)
  Function returns `false' if and only if objFunc->evaluate() returns `false'.
*/
bool CompositeOFScalarMultiply::concrete_evaluate(PatchData &patch,
                                                  double &fval, MsqError &err){
    //if invalid return false without calculating fval
  if(!objFunc->evaluate(patch, fval, err)){
    fval = 0.0;
    return false;
    
  }
  fval*=alpha;
  return true;
}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::get_quality_metric_list()"
//!Returns the QualityMetric list assossiated with objFunc.
std::list<QualityMetric*> CompositeOFScalarMultiply::get_quality_metric_list(){
  return objFunc->get_quality_metric_list();
}
	
