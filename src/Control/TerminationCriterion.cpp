// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
// DESCRIPTION:
// ============
/*! \file TerminationCriterion.cpp
  
    \brief  Member functions of the Mesquite::TerminationCriterion class

    \author Michael Brewer
    \date   Feb. 14, 2003
 */

#include "TerminationCriterion.hpp"
using namespace Mesquite;
//temporary include
#include <iostream.h>
#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::TerminationCriterion"

TerminationCriterion::TerminationCriterion() 
{
  terminationCriterionFlag=NONE;
  cullingMethodFlag=NONE;
  totalFlag=NONE;
  cullingEps=0.0;
  initialOFValue=0.0;
  previousOFValue=0.0;
  lowerOFBound=0.0;
  initialGradNorm=0.0;
    //itial size of the gradient array
  gradSize=10;
  gradNormAbsoluteEps=0.0;
  gradNormRelativeEps=0.0;
  qualityImprovementAbsoluteEps=0.0;
  qualityImprovementRelativeEps=0.0;
  iterationBound=0;
  iterationCounter=0;
  timeBound=0.0;
  vertexMovementAbsoluteEps=0.0;
  vertexMovementRelativeEps=0.0;
  successiveImprovementsAbsoluteEps=0.0;
  successiveImprovementsRelativeEps=0.0;
  boundedVertexMovementEps=0.0;
}



#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::add_criterion_type_with_double"
/*! Function to add a type of termination criterion to this object.  It is
  only valid if the specified criterion type requires a single double
  value.*/
void TerminationCriterion::add_criterion_type_with_double(TCType tc_type,
                                                       double eps,
                                                       MsqError &err)
{
  switch(tc_type){
    case GRADIENT_NORM_ABSOLUTE:
       terminationCriterionFlag|=GRADIENT_NORM_ABSOLUTE;
       gradNormAbsoluteEps=eps;
       break;
    case GRADIENT_NORM_RELATIVE:
       terminationCriterionFlag|=GRADIENT_NORM_RELATIVE;
       gradNormRelativeEps=eps;
       break;  
    case QUALITY_IMPROVEMENT_ABSOLUTE:
       terminationCriterionFlag|=QUALITY_IMPROVEMENT_ABSOLUTE;
       qualityImprovementAbsoluteEps=eps;
       break;
    case QUALITY_IMPROVEMENT_RELATIVE:
       terminationCriterionFlag|=QUALITY_IMPROVEMENT_RELATIVE;
       qualityImprovementRelativeEps=eps;
       break;
    case CPU_TIME:
       terminationCriterionFlag|=CPU_TIME;
       timeBound=eps;
       break;
    case VERTEX_MOVEMENT_ABSOLUTE:
       terminationCriterionFlag|=VERTEX_MOVEMENT_ABSOLUTE;
         //we actually compare squared movement to squared epsilon
       vertexMovementAbsoluteEps=(eps*eps);
       break;
    case VERTEX_MOVEMENT_RELATIVE:
       terminationCriterionFlag|=VERTEX_MOVEMENT_RELATIVE;
         //we actually compare squared movement to squared epsilon
       vertexMovementRelativeEps=(eps*eps);
       break;
    case SUCCESSIVE_IMPROVEMENTS_ABSOLUTE:
       terminationCriterionFlag|=SUCCESSIVE_IMPROVEMENTS_ABSOLUTE;
       successiveImprovementsAbsoluteEps=eps;
       break;
    case SUCCESSIVE_IMPROVEMENTS_RELATIVE:
       terminationCriterionFlag|=SUCCESSIVE_IMPROVEMENTS_RELATIVE;
       successiveImprovementsRelativeEps=eps;
       break; 
    default:
       err.set_msg("TCType not valid for this function.");
  };
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::add_criterion_type_with_int"
/*! Function to add a type of termination criterion to this object.  It is
  only valid if the specified criterion type requires a single integer
  value.*/
void TerminationCriterion::add_criterion_type_with_int(TCType tc_type,
                                                    int bound,
                                                    MsqError &err)
{
  switch(tc_type){
    case NUMBER_OF_ITERATES:
       terminationCriterionFlag|=NUMBER_OF_ITERATES;
       iterationBound=bound;
       break;
    default:
       err.set_msg("TCType not valid for this function.");
  };
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::set_culling_type"
/*! Function to add a type of termination criterion to this object which
  will be used for culling purposes only.  As part of the global
  termination criterion, the outer criterion will also check to make
  sure that there are still free vertices in the mesh.*/
void TerminationCriterion::set_culling_type(TCType tc_type, double eps,
                                         MsqError &err)
{
  switch(tc_type){
    case QUALITY_IMPROVEMENT_ABSOLUTE:
       cullingMethodFlag=QUALITY_IMPROVEMENT_ABSOLUTE;
       break;
    case QUALITY_IMPROVEMENT_RELATIVE:
       cullingMethodFlag=QUALITY_IMPROVEMENT_RELATIVE;
       break;
    case VERTEX_MOVEMENT_ABSOLUTE:
       cullingMethodFlag=VERTEX_MOVEMENT_ABSOLUTE;
       break;
    case VERTEX_MOVEMENT_RELATIVE:
       cullingMethodFlag=VERTEX_MOVEMENT_RELATIVE;
       break;   
    case SUCCESSIVE_IMPROVEMENTS_ABSOLUTE:
       cullingMethodFlag=SUCCESSIVE_IMPROVEMENTS_ABSOLUTE;
       break;
     case SUCCESSIVE_IMPROVEMENTS_RELATIVE:
       cullingMethodFlag=SUCCESSIVE_IMPROVEMENTS_RELATIVE;
       break;  
    default:
       err.set_msg("TCType not valid for this function.");
  };
  cullingEps=eps;
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::initialize"
/*!
  PatchData &pd is only in the arguement list so that we can
  reset the memento (???).  So, it should be the PatchData
  object that is going to be used in the optimization.
  TODO:  When culling, we  need to reset the soft_fixed flag to
  free either here or in cleanup();
 */
void TerminationCriterion::initialize(MeshSet &/*ms*/, PatchData &pd,
                                      MsqError &err)
{
    //compute total bit flag
  totalFlag = (terminationCriterionFlag | cullingMethodFlag);
  
  if(totalFlag & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE) ){
      //we need to reset the vertex memento. Do we need to keep
      //in mind that the Termination Criterion Object may have been used
      //previously?  We must delete the memento that we create here.
    previousVerticesMemento=pd.create_vertices_memento(err);
    if(terminationCriterionFlag & (VERTEX_MOVEMENT_RELATIVE)){
      initialVerticesMemento=pd.create_vertices_memento(err);
    }
    MSQ_CHKERR(err);
  }
    //if needed create an array so that it will be there
  if(totalFlag & ( GRADIENT_NORM_ABSOLUTE | GRADIENT_NORM_RELATIVE ) ){
    mGrad = new Vector3D[gradSize];
  }
    
  
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::reset"
/*!This version of reset is called using a MeshSet, which implies
  it is only called when this criterion is used as the 'outer' termination
  criterion.  
 */
void TerminationCriterion::reset(MeshSet &ms, ObjectiveFunction* obj_ptr,
                                 MsqError &err)
{
  PatchData global_patch;
    //if we need to fill out the global patch data object.
  if(totalFlag & (QUALITY_IMPROVEMENT_ABSOLUTE | QUALITY_IMPROVEMENT_RELATIVE
                  | SUCCESSIVE_IMPROVEMENTS_ABSOLUTE
                  | SUCCESSIVE_IMPROVEMENTS_RELATIVE
                  | GRADIENT_NORM_ABSOLUTE
                  | GRADIENT_NORM_RELATIVE  ) ){
    globalPatchParams.set_patch_type(PatchData::GLOBAL_PATCH, err,0,0);
    ms.get_next_patch(global_patch,globalPatchParams,err);
  }
    //currently set an error if the user is trying to terminate on the
    //outer loop using vertex movement
  if((totalFlag) & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE)){
    err.set_msg("Outer loop termination criterion on vertex movement, not implemented.");
  }
    //now call the other reset
  reset(global_patch,obj_ptr,err);MSQ_CHKERR(err);
}

    
#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::reset"
/*!

 */
void TerminationCriterion::reset(PatchData &pd, ObjectiveFunction* obj_ptr,
                                    MsqError &err)
{
    //reset the inner counter if needed
  if(totalFlag & NUMBER_OF_ITERATES){
    iterationCounter=0;
  }
    //recreate the initial memento if needed
  if(totalFlag & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE) ){
      //we need to store the previous vertex memento.
    pd.recreate_vertices_memento(previousVerticesMemento,err);
    if(totalFlag & (VERTEX_MOVEMENT_RELATIVE)){
        //we need to store the initial vertex memento.
      pd.recreate_vertices_memento(initialVerticesMemento,err);
    }
  }
    //reset the inner timer if needed
  if(totalFlag & CPU_TIME){
    mTimer.reset();
  }
    //find the initial objective function value if needed
  if(totalFlag & (QUALITY_IMPROVEMENT_ABSOLUTE | QUALITY_IMPROVEMENT_RELATIVE
                  | SUCCESSIVE_IMPROVEMENTS_ABSOLUTE
                  | SUCCESSIVE_IMPROVEMENTS_RELATIVE ) ){
    if(!obj_ptr->evaluate(pd, initialOFValue, err)){
      err.set_msg("Initial patch is invalid for evaluation.");
      MSQ_CHKERR(err);
    }
      //std::cout<<"\nReseting initial of value = "<<initialOFValue;
    previousOFValue=initialOFValue;
  }
    //GRADIENT
  if(totalFlag & (GRADIENT_NORM_ABSOLUTE | GRADIENT_NORM_RELATIVE))
  {
    int num_vertices=pd.num_vertices();
      //if the array, mGrad, is not large enough, lengthen it
    if(num_vertices>gradSize){
      delete []mGrad;
      gradSize=num_vertices;
      mGrad = new Vector3D[gradSize];
    }
      //if the initial gradient norm is needed
    if(totalFlag & (GRADIENT_NORM_RELATIVE)){
        //get gradient and make sure it is valid
      if(!obj_ptr->compute_gradient(pd, mGrad , err, num_vertices)){
        err.set_msg("Initial patch is invalid for gradient compuation.");
      }
        //get the norm of the gradient
      initialGradNorm = compute_gradient_norm(mGrad, num_vertices, err);
        //std::cout<<"\nInitial gradient norm = "<<initialGradNorm;
    }
  }   
  
}


#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::terminate"
/*!  This function evaluates the needed information and then evaluates
  the termination criteria.  If any of the selected criteria are satisfied,
  the function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::terminate(PatchData &pd, ObjectiveFunction* obj_ptr,
                                   MsqError &err)
{
    //std::cout<<"\nInside terminate(pd,of,err):  flag = "<<terminationCriterionFlag;
  
    //if terminating on numbering of inner iterations
  if(terminationCriterionFlag & NUMBER_OF_ITERATES){
    ++iterationCounter;
      //std::cout<<"\nInside terminate(pd,of,err):  "<<iterationCounter<<"  less than "<<iterationBound;
    
    if(iterationCounter>=iterationBound){
      return true;
    }
  }
    //if terminating on inner cpu time
  if(terminationCriterionFlag & CPU_TIME){
      //std::cout<<"\nCHECKED CPU time value = "<<mTimer.since_birth();
    if(mTimer.since_birth()>=timeBound){
      return true;
    }
  }
    //if terminating on vertex movement
  if(terminationCriterionFlag & ( VERTEX_MOVEMENT_ABSOLUTE |
                                  VERTEX_MOVEMENT_ABSOLUTE) ){
      //we need a PatchData memeber funciton which returns the
      //max distance (squared) of current vertex positions to the
      //old
    double max_movement_sqr=
       pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
      //std::cout<<"\nNODE_MOVEMENT = "<<max_movement_sqr;
    if(terminationCriterionFlag & ( VERTEX_MOVEMENT_ABSOLUTE) ){
      if(max_movement_sqr<=vertexMovementAbsoluteEps){
        return true;
      }
    }
    if(terminationCriterionFlag & ( VERTEX_MOVEMENT_RELATIVE) ){
      double max_movement_sqr_abs=
         pd.get_max_vertex_movement_squared(initialVerticesMemento,err);
      if(max_movement_sqr <= (vertexMovementRelativeEps*max_movement_sqr_abs))
      {
        return true;
      }
    } 
      //else we need to store the new memento
    pd.recreate_vertices_memento(previousVerticesMemento,err);
  }
    //if we need to compute the objective function value
  if(terminationCriterionFlag & (QUALITY_IMPROVEMENT_ABSOLUTE |
                                 QUALITY_IMPROVEMENT_RELATIVE |
                                 SUCCESSIVE_IMPROVEMENTS_ABSOLUTE |
                                 SUCCESSIVE_IMPROVEMENTS_RELATIVE) ){
    double obj_val;
    if(!obj_ptr->evaluate(pd, obj_val, err)){
      err.set_msg("Invalid patch passed to TerminationCriterion.");
      MSQ_CHKERR(err);
    }
      //std::cout<<"\nOF current "<<obj_val;
      //std::cout<<"\nOF previous - current "<<previousOFValue-obj_val;
      // if termination on quality improvement absolute
    if(terminationCriterionFlag & QUALITY_IMPROVEMENT_ABSOLUTE){
        //if the improvement was enough
      if(obj_val <= qualityImprovementAbsoluteEps){
        return true;
      }
    }
     // if termination on quality improvement relative
    if(terminationCriterionFlag & QUALITY_IMPROVEMENT_RELATIVE){
        //if the improvement was enough
      if((obj_val-lowerOFBound)<=(qualityImprovementRelativeEps*(
         initialOFValue-lowerOFBound))){
        return true;
      }
    }
      //if termination on successive improvements absolute
    if(terminationCriterionFlag & SUCCESSIVE_IMPROVEMENTS_ABSOLUTE){
        //if the last improvement was not significant enough
      if((previousOFValue-obj_val) <= successiveImprovementsAbsoluteEps){
        return true;
      }
    }
    //if termination on successive improvements relative
    if(terminationCriterionFlag & SUCCESSIVE_IMPROVEMENTS_RELATIVE){
        //if the last improvement was not significant enough
      if((previousOFValue-obj_val)<=(successiveImprovementsRelativeEps*(
         initialOFValue-obj_val))){
        return true;
      }
    }
      //if not termination, update the previousOFValue to by obj_val
    previousOFValue=obj_val;
  }//end of termination criteria which need objective calculated

    //if terminating on the norm of the gradient
  if(terminationCriterionFlag & (GRADIENT_NORM_ABSOLUTE |
                                 GRADIENT_NORM_RELATIVE) ){
    int num_vertices=pd.num_vertices();
      //if the array, mGrad, is too small, set an error for now
    if(num_vertices>gradSize){
      err.set_msg("Number of vertices has increased making the gradient too small");
    }
      //get gradient and make sure it is valid
    if(!obj_ptr->compute_gradient(pd, mGrad , err, num_vertices)){
      err.set_msg("Initial patch is invalid for gradient compuation.");
    }
      //get the norm of the gradient
    double grad_norm = compute_gradient_norm(mGrad, num_vertices, err);
      //std::cout<<"\n Gradient norm = "<<grad_norm;
      //if stopping on norm of the gradient absolute
    if(terminationCriterionFlag & GRADIENT_NORM_ABSOLUTE){
      if(grad_norm <= gradNormAbsoluteEps){
        return true;
      }
    }
      //if stopping on norm of the gradient relative
    if(terminationCriterionFlag & GRADIENT_NORM_RELATIVE){
      if(grad_norm <= (gradNormRelativeEps*initialGradNorm))
      {
          //std::cout<<"\nALKSJF eps ="<<gradNormRelativeEps;
          //std::cout<<"\nALKSJ init ="<<initialGradNorm;;
        return true;
      }
    }
  }
    //if none of the criteria were satisfied
  return false;
}
  

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::terminate"
/*!  This version of terminate only creates a global patch from the
  given MeshSet and then calls terminate using that PatchData object.
  Currently, this function sets an error if VERTEX_MOVEMENT is one of
  the selected criteria, because that case has not yet been implemented.

 */
bool TerminationCriterion::terminate(MeshSet &ms, ObjectiveFunction* obj_ptr,
                             MsqError &err)
{
  PatchData global_patch;
    //if we need to fill out the global patch data object.
  if(terminationCriterionFlag & (QUALITY_IMPROVEMENT_ABSOLUTE |
                                 QUALITY_IMPROVEMENT_RELATIVE |
                                 SUCCESSIVE_IMPROVEMENTS_ABSOLUTE |
                                 SUCCESSIVE_IMPROVEMENTS_RELATIVE |
                                 GRADIENT_NORM_ABSOLUTE |
                                 GRADIENT_NORM_RELATIVE )){
    globalPatchParams.set_patch_type(PatchData::GLOBAL_PATCH, err,0,0);
    
    ms.get_next_patch(global_patch,globalPatchParams,err);
  }
    //currently set an error if the user is trying to terminate on the
    //outer loop using vertex movement
  if(terminationCriterionFlag & (VERTEX_MOVEMENT_ABSOLUTE |
                                 VERTEX_MOVEMENT_RELATIVE ) ){
    err.set_msg("Outer loop termination criterion on vertex movement, not implemented.");
  }
    //now call the other terminate... with the patchdata object.
  bool return_bool=terminate(global_patch,obj_ptr,err);
  MSQ_CHKERR(err);
  return return_bool;
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::cull_vertices"
/*!This function check the culling method criterion supplied to the object
  by the user.  If the user does not supply a culling method criterion,
  the default criterion is NONE, and in that case, no culling is performed.
  If the culling method criterion is satisfied, the interior vertices
  of the given patch are flagged as soft_fixed.  Otherwise, the soft_fixed
  flag is removed from each of the vertices in the patch (interior or
  boundary vertices).  Also, if the criterion was satisfied, then the
  function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::cull_vertices(PatchData &pd,
                                      ObjectiveFunction* obj_ptr,
                                      MsqError &err)
{
  bool cull_bool=false;
  double obj_val=0.0;
  switch(cullingMethodFlag){
    case NONE:
       return cull_bool;
    case QUALITY_IMPROVEMENT_ABSOLUTE:
       if(!obj_ptr->evaluate(pd, obj_val, err)){
         err.set_msg("Invalid patch passed to TerminationCriterion.");
         MSQ_CHKERR(err);
       }
         //if the improvement was enough
       if((obj_val)<=(cullingEps))
       {
         cull_bool=true;  
       }
       else
       {
         cull_bool=false;
       }
       break;
    case QUALITY_IMPROVEMENT_RELATIVE:
       if(!obj_ptr->evaluate(pd, obj_val, err)){
         err.set_msg("Invalid patch passed to TerminationCriterion.");
         MSQ_CHKERR(err);
       }
         //if the improvement was enough
       if((obj_val-lowerOFBound)<=(cullingEps*(initialOFValue-lowerOFBound)))
       {
         cull_bool=true;  
       }
       else
       {
         cull_bool=false;
       }
       break;
    case VERTEX_MOVEMENT_ABSOLUTE:
         //if movement was small enough
       if(pd.get_max_vertex_movement_squared(previousVerticesMemento,err)<
          cullingEps){
         cull_bool=true;  
       }
       else
       {
         cull_bool=false;
       }
       break;
    case VERTEX_MOVEMENT_RELATIVE:
         //if movement was small enough
       if(pd.get_max_vertex_movement_squared(previousVerticesMemento,err)<
          (cullingEps*
           pd.get_max_vertex_movement_squared(initialVerticesMemento,err))){
         cull_bool=true;  
       }
       else
       {
         cull_bool=false;
       }
       break;
    default:
       err.set_msg("Requested culling method not yet implemented.");
  };
    //Now actually have patch data cull vertices
  if(cull_bool){
    pd.set_all_vertices_soft_free(err);
  }
  else{
    pd.set_free_vertices_soft_fixed(err);
  }
  MSQ_CHKERR(err);
  return cull_bool;
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::cleanup"
/*!
  Currently this only deletes the memento of the vertex positions and the
  mGrad vector if neccessary.
  TODO:  When culling, we  need to reset the soft_fixed flag to
  free either here or in initialize();
 */
void TerminationCriterion::cleanup(MeshSet &/*ms*/, MsqError &/*err*/)
{
    //clean up the memento if used
  if(totalFlag & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE )){
    previousVerticesMemento->~PatchDataVerticesMemento();
    if(terminationCriterionFlag & (VERTEX_MOVEMENT_RELATIVE))
    {
      initialVerticesMemento->~PatchDataVerticesMemento();
    }
  }
  if(totalFlag & (GRADIENT_NORM_ABSOLUTE | GRADIENT_NORM_RELATIVE ) ){
    delete [] mGrad;
  } 
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::cleanup"
/*!
  Currently this function only computes the infinity norm of the
  gradient array, vec.  .
  TODO:  This function needs to be able to compute other norms.
 */
double TerminationCriterion::compute_gradient_norm(Vector3D * const vec,
                                                   int len, MsqError &/*err*/)
{
  double grad_norm=0;
  for(int gi=0;gi<len;++gi){
    for (int gj=0;gj<3;++gj){
      if(grad_norm<fabs(vec[gi][gj])){
        grad_norm=fabs(vec[gi][gj]);
      }
    }
  }
  return grad_norm;
}
