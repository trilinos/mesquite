/*!
  \file    CompositeOFMultiply.cpp
  \brief  

  This Objective Function combines two Objective Functions by mulitplication
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "ObjectiveFunction.hpp"
#include "CompositeOFMultiply.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "CompositeOFMultiply::CompositeOFMultiply"
/*!
Sets the QualityMetric pointer to the metric associated with Obj1 and Obj2
if Obj1 and Obj2 are associated with the same metric.  Otherwise, it sets
the QualityMetric pointer to NULL.  If Obj1 or Obj2 requires a feasible
region, then so does the new CompositeOFMultiply ObjectiveFunction.  The new
ObjectiveFunction's negateFlag is set to negative one only if both Obj1 and
Obj2's negateFlag are negative one (because obj1 and obj2's evaluate function
multiply their return values by negative one if their respective
function needs to be maximized.  If both of these functions needed to
be maximized, then the negative ones will have cancelled out).  Otherwise,
the negateFlag is set to one.
  \param Obj1 (ObjectiveFunction*)
  \param Obj2 (ObjectiveFunction*)
 */
CompositeOFMultiply::CompositeOFMultiply(ObjectiveFunction* Obj1, ObjectiveFunction* Obj2){
  if(Obj1->get_quality_metric()==Obj2->get_quality_metric()){
    set_quality_metric(Obj1->get_quality_metric());
  }
  else
    set_quality_metric(NULL);
  if(Obj1->get_feasible_constraint())
    set_feasible(Obj1->get_feasible_constraint());
  else
    set_feasible(Obj2->get_feasible_constraint());
  objFunc1=Obj1;
  objFunc2=Obj2;
    //if both obj1 and ob2 have been negated
  if(Obj1->get_negate_flag()-Obj2->get_negate_flag()==-2)
    set_negate_flag(-1);
  else
    set_negate_flag(1);
}

#undef __FUNC__
#define __FUNC__ "CompositeOFMultiply::~CompositeOFMultiply"

//Michael:  need to clean up here
CompositeOFMultiply::~CompositeOFMultiply(){

}

#undef __FUNC__
#define __FUNC__ "CompositeOFMultiply::get_quality_metric_list"
/*! Returns the QualityMetric list associated with objFunc1 merged
  with the QualityMetric list associated with objFunc2.  The entries
  in this merged list may not be unique.
*/
std::list<QualityMetric*> CompositeOFMultiply::get_quality_metric_list()
{
  std::list<QualityMetric*> temp_list=objFunc1->get_quality_metric_list();
  std::list<QualityMetric*> temp_list2=objFunc2->get_quality_metric_list();
  temp_list.merge(temp_list2);
  return temp_list;
}


#undef __FUNC__
#define __FUNC__ "CompositeOFMultiply::concrete_evaluate"
/*!Returns objFunc1->evaluate(patch,err)*objFunc2->evaluate(patch,err).
  Note that since objFunc1 and objFunc2's evaluate() functions are called
  (as opposed to their concrete_evaluates) the returned values have
  already been multiplied by the respective negateFlag (that is,
  if objFunc1 (or objFunc2) needed to be maximized then the value has
  been multiplied by negative one so that it may be minimized instead.)
*/
double CompositeOFMultiply::concrete_evaluate(PatchData &patch, MsqError &err){
  //Total value of objective function
  double total_value=objFunc1->evaluate(patch, err);
  total_value*=objFunc2->evaluate(patch, err);
  return total_value;
}
	
	
