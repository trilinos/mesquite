/*!
  \file    CompositeOFAdd.cpp
  \brief  

  This Objective Function combines two Objective Functions by addition
  \author Michael Brewer
  \date   2002-06-24
*/
#include <math.h>
#include "ObjectiveFunction.hpp"
#include "CompositeOFAdd.hpp"
#include "MsqTimer.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "CompositeOFAdd::CompositeOFAdd"
/*!
Sets the QualityMetric pointer to the metric associated with Obj1 and Obj2
if Obj1 and Obj2 are associated with the same metric.  Otherwise, it sets
the QualityMetric pointer to NULL.  The new
ObjectiveFunction's negateFlag is always set to one, because the values
produced by obj1 and obj2 have already been multiplied by negative one
if it was needed.  Defaults to the analytical gradient.
  \param Obj1 (ObjectiveFunction*)
  \param Obj2 (ObjectiveFunction*)
 */
CompositeOFAdd::CompositeOFAdd(ObjectiveFunction* Obj1,
                               ObjectiveFunction* Obj2){
  if(Obj1->get_quality_metric()==Obj2->get_quality_metric()){
    set_quality_metric(Obj1->get_quality_metric());
  }
  else
    set_quality_metric(NULL);
  objFunc1=Obj1;
  objFunc2=Obj2;
  set_negate_flag(1);
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
}

#undef __FUNC__
#define __FUNC__ "CompositeOFAdd::~CompositeOFAdd"

//Michael:  need to clean up here
CompositeOFAdd::~CompositeOFAdd(){

}

#undef __FUNC__
#define __FUNC__ "CompositeOFAdd::get_quality_metric_list"
/*! Returns the QualityMetric list associated with objFunc1 merged
  with the QualityMetric list associated with objFunc2.  The entries
  in this merged list may not be unique.
*/
list<QualityMetric*> CompositeOFAdd::get_quality_metric_list()
{
  list<QualityMetric*> temp_list=objFunc1->get_quality_metric_list();
  list<QualityMetric*> temp_list2=objFunc2->get_quality_metric_list();
  temp_list.merge(temp_list2);
  return temp_list;
    
}

#undef __FUNC__
#define __FUNC__ "CompositeOFAdd::concrete_evaluate"
/*!Compute fval= objFunc1->evaluate(patch,err)+objFunc2->evaluate(patch,err).
  Note that since objFunc1 and objFunc2's evaluate() functions are called
  (as opposed to their concrete_evaluates) the returned values have
  already been multiplied by the respective negateFlag (that is,
  if objFunc1 (or objFunc2) needed to be maximized then the value has
  been multiplied by negative one so that it may be minimized instead.)
  Function returns `false' if  either objFunc1->concrete_evaluate() or
  objFunc2->concrete_evaluate() returns `false'; otherwise, function
  returns `true'.
*/
bool CompositeOFAdd::concrete_evaluate(PatchData &patch, double &fval,
                                       MsqError &err){
  double second_val=0.0;
    //If patch is invalid
  if(! objFunc1->evaluate(patch, fval, err)){
    fval=0.0;
    return false;
  }
  if(! objFunc2->evaluate(patch, second_val, err)){
    fval=0.0;
    return false;
  }
  fval+=second_val;
  return true;
}
	
	
#undef __FUNC__
#define __FUNC__ "CompositeOFAdd::compute_analytical_gradient"
/*! Analytically computes the composite objective function's gradient
  by combining the gradients returned from 
  objFunc2->compute_gradient() and objFunc2->compute_gradient().
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param OF_val The objective function value.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool CompositeOFAdd::compute_analytical_gradient(PatchData &patch,
                                                 Vector3D *const &grad,
						 double &OF_val,
                                                 MsqError &err,
                                                 size_t array_size)
{
  FUNCTION_TIMER_START(__FUNC__);
  double second_val=0.0;//store the second objective function val
  OF_val=0.0;
    //get first objective function's gradient
  bool rval=objFunc1->compute_gradient(patch, grad, OF_val, 
				       err, array_size);
  if(rval){
    int num_vert=patch.num_vertices();
    Vector3D* second_grad = new Vector3D[num_vert];
      //get second objective function's gradient
    rval=objFunc2->compute_gradient(patch, second_grad, second_val,
				    err, num_vert);
      //if both objective functions were successfully computed, add them
    if(rval){
      //add the two objective function values.
      OF_val+=second_val;
      int i=0;
      for(i=0;i<num_vert;++i){
        grad[i]+=second_grad[i];
      }
        //delete the dynamically allocated space for the second gradient
      delete []second_grad;
    }
  }
  FUNCTION_TIMER_END();
    //true if both of the above compute gradient's were successful.
  if(!rval)
    OF_val=0.0;
  return rval;
}
