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
#include "MsqTimer.hpp"
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
the negateFlag is set to one.  Defaults to the analytical gradient.
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
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
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
/*!Computes fval= objFunc1->evaluate(patch,err)*objFunc2->evaluate(patch,err).
  Note that since objFunc1 and objFunc2's evaluate() functions are called
  (as opposed to their concrete_evaluates) the returned values have
  already been multiplied by the respective negateFlag (that is,
  if objFunc1 (or objFunc2) needed to be maximized then the value has
  been multiplied by negative one so that it may be minimized instead.)
  Function returns `false' if either objFunc1->evaluate() or
  objFunc2->evaluate() returns `false'; otherwise function returns `true'.
*/
bool CompositeOFMultiply::concrete_evaluate(PatchData &patch, double &fval,
                                            MsqError &err){
  double second_val;
    //if invalid, return false without calculating fval.
  if(! objFunc1->evaluate(patch, fval, err)){
    fval=0.0;
    return false;
  }
  
  if(! objFunc2->evaluate(patch, second_val, err)){
    fval=0.0;
    return false;
  }
  fval*=second_val;
  return true;
}

#undef __FUNC__
#define __FUNC__ "CompositeOFMultiply::compute_analytical_gradient"
/*! Analytically computes the composite objective function's gradient,
  using the multiplication rule.
  by scaling the gradient returned objFunc->compute_gradient().
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool CompositeOFMultiply::compute_analytical_gradient(PatchData &patch,
                                                 Vector3D *const &grad,
                                                 MsqError &err,
                                                 int array_size)
{
#ifdef USE_FUNCTION_TIMERS          
  StopWatchCollection::Key this_key = GlobalStopWatches.add(__FUNC__,false);
  GlobalStopWatches.start(this_key);
#endif
  double obj_1_val=0.0;
  double obj_2_val=0.0;
    //get first obj function's value
  bool rval=objFunc1->evaluate(patch, obj_1_val, err);
  if(rval){
      //if the first value was valid, get the second
    rval=objFunc2->evaluate(patch, obj_2_val, err);
  }
    
  if(rval){
      //if both obj function values were valid, get first gradient
    rval=objFunc1->compute_gradient(patch, grad, err, array_size);
      //if all of the above were valid, get the second gradient
    if(rval){
      int num_vert=patch.num_vertices();
      Vector3D* second_grad = new Vector3D[num_vert];
        //get second objective function's gradient
      rval=objFunc2->compute_gradient(patch, second_grad, err, num_vert);
        //if both objective functions were successfully computed, and
        //both evaluates were successful, use the multiplaction rule
        //to get the complete gradient.
      if(rval){
        int i=0;
        for(i=0;i<num_vert;++i){
          grad[i]*=obj_2_val;
          grad[i]+=(second_grad[i]*obj_1_val);
        }
      }
        //delete the dynamically allocated space for the second gradient
      delete []second_grad;
    }
  }
#ifdef USE_FUNCTION_TIMERS          
  GlobalStopWatches.stop(this_key);
#endif
    //true if both gradient and both evaluate were successful.
  return rval;
}
	
	
