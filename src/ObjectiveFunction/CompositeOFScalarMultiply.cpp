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
ObjectiveFunction's negateFlag is the opposite of Obj's.  This objective
function defaults to the analytical gradient.
  \param alp (double)
  \param Obj (ObjectiveFunction*)
 */
CompositeOFScalarMultiply::CompositeOFScalarMultiply(double alp, ObjectiveFunction* Obj){

  set_quality_metric(Obj->get_quality_metric());
  set_feasible(Obj->get_feasible_constraint());
  objFunc=Obj;
  mAlpha=alp;
  if(alp<0)
    set_negate_flag(-1);
  else if (alp>0)
    set_negate_flag(1);
  else
    PRINT_WARNING("ObjectiveFunction being scaled by zero.");
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::~CompositeOFScalarMultiply"

//Michael:  need to clean up here
CompositeOFScalarMultiply::~CompositeOFScalarMultiply(){

}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::concrete_evaluate"
/*!Computes fval= mAlpha*objFunc->evaluate(patch,err).  Note that since Obj's
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
  fval*=mAlpha;
  return true;
}

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::get_quality_metric_list()"
//!Returns the QualityMetric list assossiated with objFunc.
std::list<QualityMetric*> CompositeOFScalarMultiply::get_quality_metric_list(){
  return objFunc->get_quality_metric_list();
}
	

#undef __FUNC__
#define __FUNC__ "CompositeOFScalarMultiply::compute_analytical_gradient"
/*! Analytically computes the composite objective function's gradient
  by scaling the gradient returned objFunc->compute_gradient().  If
  mAlpha is less than zero, the gradient is scaled by negatvie mAlpha
  because the ObjectiveFunction is multiplied by negative one so that
  it may be minimized.
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool CompositeOFScalarMultiply::compute_analytical_gradient(PatchData &patch,
                                                            Vector3D *const
                                                            &grad,
                                                            MsqError &err,
                                                            int array_size)
{
#ifdef USE_FUNCTION_TIMERS          
  StopWatchCollection::Key this_key = GlobalStopWatches.add(__FUNC__,false);
  GlobalStopWatches.start(this_key);
#endif
  double scale_factor=(get_negate_flag()*mAlpha);
  bool rval=objFunc->compute_gradient(patch, grad, err, array_size);
  int num_vert=patch.num_vertices();
  int i=0;
    //If the objFunc was successful in calculating the gradient
  if(rval){
      //scale the gradient by alpha
    for(i=0;i<num_vert;++i){
      grad[i]*=scale_factor;
    }
  }
#ifdef USE_FUNCTION_TIMERS          
  GlobalStopWatches.stop(this_key);
#endif
  return rval;
}
