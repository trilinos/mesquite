/*!
  \file    CompositeOFAdd.cpp
  \brief  

  This Objective Function combines two Objective Functions by addition
  \author Michael Brewer
  \date   2002-06-24
*/
#include <math.h>
#include "ObjectiveFunction.hpp"
#include "CompositeOFScalarAdd.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::CompositeOFScalarAdd"
/*!
Sets the QualityMetric pointer to the metric associated with Obj1.  If
Obj1 requires a feasible region, then so does the new CompositeOFScalarAdd
ObjectiveFunction.  The new ObjectiveFunction's negateFlag is also the
same as that of Obj1..  
  \param alp (double)
  \param Obj1 (ObjectiveFunction*)
 */
CompositeOFScalarAdd::CompositeOFScalarAdd(double alp, ObjectiveFunction* Obj1){
  set_quality_metric(Obj1->get_quality_metric());
  set_feasible(Obj1->get_feasible_constraint());
  objFunc=Obj1;
  alpha=alp;
  set_negate_flag(1);
}


#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::~CompositeOFScalarAdd"

//Michael:  need to clean up here
CompositeOFScalarAdd::~CompositeOFScalarAdd(){

}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::concrete_evaluate()"
/*!Returns alpha+objFunc->evaluate(patch,err).  Note that since Obj's
  evaluate() function is called (as opposed to its concrete_evaluate) the
  returned value has been multiplied by objFunc's negateFlag (that is,
  if objFunc needed to be maximized then the value has been multiplied
  by negative one so that it may be minimized instead.)
*/
double CompositeOFScalarAdd::concrete_evaluate(PatchData &patch, MsqError &err){
  //Total value of objective function
  double total_value=objFunc->evaluate(patch, err);
    //if(objFunc2)
    //total_value*=objFunc2->evaluate(patch, err);
  total_value+=alpha;
  return total_value;
}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::get_quality_metric_list()"
//!Returns the QualityMetric list assossiated with objFunc.
std::list<QualityMetric*> CompositeOFScalarAdd::get_quality_metric_list(){
  return objFunc->get_quality_metric_list();
}


