/*!
  \file    CompositeOFScalarAdd.cpp
  \brief  

  This Objective Function combines two Objective Functions by addition
  \author Michael Brewer
  \date   2002-06-24
*/
#include <math.h>
#include "ObjectiveFunction.hpp"
#include "CompositeOFScalarAdd.hpp"
#include "MsqTimer.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::CompositeOFScalarAdd"
/*!
Sets the QualityMetric pointer to the metric associated with Obj1.  
The new ObjectiveFunction's negateFlag is also the
same as that of Obj1.  This objective function defaults to the analytical
gradient which essentially just calls Obj1's gradient function.
  \param alp (double)
  \param Obj1 (ObjectiveFunction*)
 */
CompositeOFScalarAdd::CompositeOFScalarAdd(double alp, ObjectiveFunction* Obj1){
  set_quality_metric(Obj1->get_quality_metric());
  objFunc=Obj1;
  mAlpha=alp;
  set_negate_flag(1);
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
}


#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::~CompositeOFScalarAdd"

//Michael:  need to clean up here
CompositeOFScalarAdd::~CompositeOFScalarAdd(){

}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::concrete_evaluate()"
/*!Computes fval= mAlpha+objFunc->evaluate(patch,err).  Note that since Obj's
  evaluate() function is called (as opposed to its concrete_evaluate) the
  returned value has been multiplied by objFunc's negateFlag (that is,
  if objFunc needed to be maximized then the value has been multiplied
  by negative one so that it may be minimized instead.)
  Functions returns `false' if and only if objFunc->evaluate() returns
  `false'.
*/
bool CompositeOFScalarAdd::concrete_evaluate(PatchData &patch, double &fval,
                                             MsqError &err){
    //if invalid return false without calculating fval.
  if( ! objFunc->evaluate(patch, fval, err)){
    fval = 0.0;
    return false;
  }
  
  fval+=mAlpha;
  return true;
}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::get_quality_metric_list()"
//!Returns the QualityMetric list assossiated with objFunc.
list<QualityMetric*> CompositeOFScalarAdd::get_quality_metric_list(){
  return objFunc->get_quality_metric_list();
}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarAdd::compute_analytical_gradient"
/*! Analytically computes the composite objective function's gradient
  by scaling the gradient returned objFunc->compute_gradient().
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool CompositeOFScalarAdd::compute_analytical_gradient(PatchData &patch,
                                                       Vector3D *const &grad,
						       double &OF_val,
                                                       MsqError &err,
                                                       size_t array_size)
{
  FUNCTION_TIMER_START(__FUNC__);
  bool rval=objFunc->compute_gradient(patch, grad, OF_val,err, array_size);
  OF_val+=mAlpha;
  FUNCTION_TIMER_END();
  return rval;
}

