/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
// DESCRIPTION:
// ============
/*! \file TerminationCriterion.cpp
  
    \brief  Member functions of the Mesquite::TerminationCriterion class

    \author Michael Brewer
    \author Thomas Leurent
    \date   Feb. 14, 2003
 */

#include "TerminationCriterion.hpp"
#include "MeshSet.hpp"
#include "MsqVertex.hpp"
using namespace Mesquite;

#include "MesquiteInterrupt.hpp"

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::TerminationCriterion"
/*!Constructor initializes all of the data members which are not
  necessarily automatically initialized in their constructors.*/
TerminationCriterion::TerminationCriterion() 
{
  terminationCriterionFlag=NONE;
  cullingMethodFlag=NONE;
  totalFlag=NONE;
  cullingEps=0.0;
  initialOFValue=0.0;
  previousOFValue=0.0;
  lowerOFBound=0.0;
  initialGradL2Norm=0.0;
  initialGradInfNorm=0.0;
    //initial size of the gradient array
  gradSize=10;
  gradL2NormAbsoluteEps=0.0;
  gradL2NormRelativeEps=0.0;
  gradInfNormAbsoluteEps=0.0;
  gradInfNormRelativeEps=0.0;
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
    //variables for supplied information
  functionSupplied=false;
  gradientSupplied=false;
  suppliedGradientArray=NULL;
  currentOFValue = 0.0;
  
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
    case GRADIENT_L2_NORM_ABSOLUTE:
       terminationCriterionFlag|=GRADIENT_L2_NORM_ABSOLUTE;
       gradL2NormAbsoluteEps=eps;
       break; 
    case GRADIENT_INF_NORM_ABSOLUTE:
       terminationCriterionFlag|=GRADIENT_INF_NORM_ABSOLUTE;
       gradInfNormAbsoluteEps=eps;
       break;
    case GRADIENT_L2_NORM_RELATIVE:
       terminationCriterionFlag|=GRADIENT_L2_NORM_RELATIVE;
       gradL2NormRelativeEps=eps;
       break;  
    case GRADIENT_INF_NORM_RELATIVE:
       terminationCriterionFlag|=GRADIENT_INF_NORM_RELATIVE;
       gradInfNormRelativeEps=eps;
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
    case BOUNDED_VERTEX_MOVEMENT:
       terminationCriterionFlag|=BOUNDED_VERTEX_MOVEMENT;
       boundedVertexMovementEps=eps;
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
#define  __FUNC__ "TerminationCriterion::remove_criterion_type"
/*! Function to remove a previously set criterion type.*/
void TerminationCriterion::remove_criterion_type(TCType tc_type,
                                                 MsqError &/*err*/)
{
  terminationCriterionFlag&=(~tc_type);
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
#define  __FUNC__ "TerminationCriterion::remove_culling"
/*!Sets the culling type to be NONE.*/
void TerminationCriterion::remove_culling(MsqError &/*err*/)
{
  cullingMethodFlag=NONE;
}

  

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::initialize"
/*!
  PatchData &pd is only in the arguement list so that we can
  reset the memento (???).  So, it should be the PatchData
  object that is going to be used in the optimization.
  
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
    if(totalFlag & (VERTEX_MOVEMENT_RELATIVE)){
      initialVerticesMemento=pd.create_vertices_memento(err);
    }
    MSQ_CHKERR(err);
  }
    //if needed create an array so that it will be there
  if(totalFlag & ( (GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_INF_NORM_ABSOLUTE)
                   | (GRADIENT_L2_NORM_RELATIVE | GRADIENT_INF_NORM_RELATIVE)  ) ){
    mGrad = new Vector3D[gradSize];
  }
    
  
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::reset"
/*!This version of reset is called using a MeshSet, which implies
  it is only called when this criterion is used as the 'outer' termination
  criterion.  
 */
bool TerminationCriterion::reset(MeshSet &ms, ObjectiveFunction* obj_ptr,
                                 MsqError &err)
{
  PatchData global_patch;
    //if we need to fill out the global patch data object.
  if(totalFlag & (QUALITY_IMPROVEMENT_ABSOLUTE | QUALITY_IMPROVEMENT_RELATIVE
                  | SUCCESSIVE_IMPROVEMENTS_ABSOLUTE
                  | SUCCESSIVE_IMPROVEMENTS_RELATIVE
                  | (GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_INF_NORM_ABSOLUTE) 
                  | (GRADIENT_L2_NORM_RELATIVE | GRADIENT_INF_NORM_RELATIVE)   ) ){
    globalPatchParams.set_patch_type(PatchData::GLOBAL_PATCH, err,0,0);
    ms.get_next_patch(global_patch,globalPatchParams,err);
  }
    //currently set an error if the user is trying to terminate on the
    //outer loop using vertex movement
  if((totalFlag) & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE)){
    err.set_msg("Outer loop termination criterion on vertex movement, not implemented.");
  }
    //now call the other reset
  return reset(global_patch,obj_ptr,err);
}

    
#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::reset"
/*!Reset function using using a PatchData object.  This function is
  called for the inner-stopping criterion directly from the
  loop over mesh function in VertexMover.  For outer criterion,
  it is called from the reset function which takes a MeshSet object.
  This function prepares the object to be used by setting the initial
  values of some of the data members.  As examples, if needed, it resets
  the cpu timer to zero, the iteration counter to zero, and the
  initial and previous objective function values to the current
  objective function value for this patch.
  The return value for this function is similar to that of terminate().
  The function returns false if the checked criteria have not been
  satisfied, and true if they have been.  reset() only checks the
  GRADIENT_INF_NORM_ABSOLUTE, GRADIENT_L2_NORM_ABSOLUTE, and the
  QUALITY_IMPROVEMENT_ABSOLUTE criteria.  Checking these criteria
  allows the QualityImprover to skip the entire optimization if
  the initial mesh satisfies the appropriate conditions.
 */
bool TerminationCriterion::reset(PatchData &pd, ObjectiveFunction* obj_ptr,
                                    MsqError &err)
{
  bool return_flag = false;
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
   
    //GRADIENT
  if(totalFlag & ((GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_INF_NORM_ABSOLUTE)
                  | (GRADIENT_L2_NORM_RELATIVE | GRADIENT_INF_NORM_RELATIVE) ))
  {
    int num_vertices=pd.num_vertices();
      //if the array, mGrad, is not large enough, lengthen it
    if(num_vertices>gradSize){
      delete []mGrad;
      gradSize=num_vertices;
      mGrad = new Vector3D[gradSize];
    }
      //get gradient and make sure it is valid
    if(!obj_ptr->compute_gradient(pd, mGrad , currentOFValue, 
                                  err, num_vertices)){
      err.set_msg("Initial patch is invalid for gradient computation.");
    }
      //get the gradient norms
    initialGradInfNorm = Linf(mGrad, num_vertices);
    initialGradL2Norm = length(mGrad, num_vertices);
      //the OFvalue comes for free, so save it
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
      //if stopping on L2 norm of the gradient
    if(terminationCriterionFlag & GRADIENT_L2_NORM_ABSOLUTE ){
      if(initialGradL2Norm <= gradL2NormAbsoluteEps){
        return_flag = true;
      }
    }
      //if stopping on Linf norm of the gradient
    if(terminationCriterionFlag & GRADIENT_INF_NORM_ABSOLUTE ){
      if(initialGradInfNorm <= gradInfNormAbsoluteEps){
        return_flag = true;
      }
    }
  }
  
  //find the initial objective function value if needed and not already
  //computed.  If we needed the gradient, we have the OF value for free.
  if((totalFlag & (QUALITY_IMPROVEMENT_ABSOLUTE | 
                   QUALITY_IMPROVEMENT_RELATIVE | 
                   SUCCESSIVE_IMPROVEMENTS_ABSOLUTE | 
                   SUCCESSIVE_IMPROVEMENTS_RELATIVE )) &&
     !(totalFlag & (GRADIENT_L2_NORM_RELATIVE | 
                    GRADIENT_INF_NORM_RELATIVE |
                    GRADIENT_L2_NORM_ABSOLUTE | 
                    GRADIENT_INF_NORM_ABSOLUTE) )){
      //ensure the obj_ptr is not null
    if(obj_ptr==NULL){
      err.set_msg("Error termination criteria set which uses objective functions, but no objective function is available.");
    }
    
    if(!obj_ptr->evaluate(pd, currentOFValue, err)){
      err.set_msg("Initial patch is invalid for evaluation.");
      MSQ_CHKERR(err);
    }
      //std::cout<<"\nReseting initial of value = "<<initialOFValue;
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
  }
  if(totalFlag & QUALITY_IMPROVEMENT_ABSOLUTE){
    if(initialOFValue <= qualityImprovementAbsoluteEps){
      return_flag = true;
    }
  }
    //return false if we have not yet returned anything
  return return_flag;
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
  bool return_flag=false;
  //  cout<<"\nInside terminate(pd,of,err):  flag = "<<terminationCriterionFlag << endl;

    //First check for an interrupt signal
  return_flag = MesquiteInterrupt::interrupt_was_signaled();
  if(return_flag){
    err.set_msg("Optimization terminated by user.");
  }
  
    //if terminating on numbering of inner iterations
  if(terminationCriterionFlag & NUMBER_OF_ITERATES){
    ++iterationCounter;
      //std::cout<<"\nInside terminate(pd,of,err):  "<<iterationCounter<<"  less than "<<iterationBound;
    
    if(iterationCounter>=iterationBound){
      return_flag = true;
    }
  }
    //if terminating on inner cpu time
  if(terminationCriterionFlag & CPU_TIME){
      //std::cout<<"\nCHECKED CPU time value = "<<mTimer.since_birth();
    if(mTimer.since_birth()>=timeBound){
      return_flag=true;
    }
  }
    //if terminating on vertex movement
  if(terminationCriterionFlag & ( VERTEX_MOVEMENT_ABSOLUTE |
                                  VERTEX_MOVEMENT_RELATIVE) ){
      //we need a PatchData memeber funciton which returns the
      //max distance (squared) of current vertex positions to the
      //old
    double max_movement_sqr=
       pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
      //std::cout<<"\nNODE_MOVEMENT = "<<max_movement_sqr;
    if(terminationCriterionFlag & ( VERTEX_MOVEMENT_ABSOLUTE) ){
      if(max_movement_sqr<=vertexMovementAbsoluteEps){
        return_flag = true;
      }
    }
      //if terminating on RELATIVE vertex movement
    if(terminationCriterionFlag & ( VERTEX_MOVEMENT_RELATIVE) ){
      double max_movement_sqr_abs=
         pd.get_max_vertex_movement_squared(initialVerticesMemento,err);
      if(max_movement_sqr <= (vertexMovementRelativeEps*max_movement_sqr_abs))
      {
        return_flag = true;
      }
    } 
      //else we need to store the new memento
    pd.recreate_vertices_memento(previousVerticesMemento,err);
  }


  //if terminating on the norm of the gradient
  if(terminationCriterionFlag & ( GRADIENT_L2_NORM_ABSOLUTE  |
                                  GRADIENT_INF_NORM_ABSOLUTE |
                                  GRADIENT_L2_NORM_RELATIVE  |
                                  GRADIENT_INF_NORM_RELATIVE ) ){
    int num_vertices=pd.num_vertices();
      //temp_grad holds the value of mGrad.  If the gradient array is
      //supplied by the QualityImprover, mGrad may be set to point to
      //another array, but we need to save a pointer to the old array,
      //because it was created using the new operator.
    Vector3D* temp_grad=mGrad;
      //if gradient is supplied use it
    if(gradientSupplied){
      mGrad=suppliedGradientArray;
    }
      //otherwise, calculate it
    else{
        //if the array, mGrad, is too small, set an error for now
      if(num_vertices>gradSize){
        err.set_msg("Number of vertices has increased making the gradient too small");
      }
        //get gradient and make sure it is valid
      if(!obj_ptr->compute_gradient(pd, mGrad, currentOFValue, err, num_vertices)){
        err.set_msg("Initial patch is invalid for gradient compuation.");
      }
    }//end else if gradient needed to be calcuated

    double grad_L2_norm=10e6;
    if (terminationCriterionFlag & (GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_L2_NORM_RELATIVE)) {
      grad_L2_norm = length(mGrad, num_vertices); // get the L2 norm
      MSQ_DEBUG_ACTION(1, {cout << "  o TermCrit -- gradient L2 norm: " << grad_L2_norm << endl;});
    }
    double grad_inf_norm=10e6;
    if (terminationCriterionFlag & (GRADIENT_INF_NORM_ABSOLUTE | GRADIENT_INF_NORM_RELATIVE)) {
      grad_inf_norm = length(mGrad, num_vertices); // get the Linf norm
    } 
    
    //if stopping on L2 norm of the gradient
    if(terminationCriterionFlag & GRADIENT_L2_NORM_ABSOLUTE ){
      if(grad_L2_norm <= gradL2NormAbsoluteEps){
        return_flag = true;
      }
    }
    //if stopping on Linf norm of the gradient
    if(terminationCriterionFlag & GRADIENT_INF_NORM_ABSOLUTE ){
            
      if(grad_inf_norm <= gradInfNormAbsoluteEps){
        return_flag = true;
      }
    }
    //if stopping on L2 norm of the gradient relative to the previous iteration
    if(terminationCriterionFlag & GRADIENT_L2_NORM_RELATIVE) {
      if(grad_L2_norm <= (gradL2NormRelativeEps*initialGradInfNorm))
      {
        return_flag = true;
      }
    }
    //if stopping on Linf norm of the gradient relative to the previous iteration
    if(terminationCriterionFlag & GRADIENT_INF_NORM_RELATIVE) {
      if(grad_inf_norm <= (gradInfNormRelativeEps*initialGradInfNorm))
      {
        return_flag = true;
      }
    }
      //reset mGrad to temp_grad so that it may be correctly deleted
    mGrad=temp_grad;
  }



  
    //if using a criterion that needs the objective function value
  if(terminationCriterionFlag & (QUALITY_IMPROVEMENT_ABSOLUTE |
                                 QUALITY_IMPROVEMENT_RELATIVE |
                                 SUCCESSIVE_IMPROVEMENTS_ABSOLUTE |
                                 SUCCESSIVE_IMPROVEMENTS_RELATIVE) ){
      //if the function val was supplied, use it
      //otherwise calcuate and set an error if invalid
    if(!functionSupplied && (!(terminationCriterionFlag &
                               ( GRADIENT_L2_NORM_ABSOLUTE  |
                                 GRADIENT_INF_NORM_ABSOLUTE |
                                 GRADIENT_L2_NORM_RELATIVE  |
                                 GRADIENT_INF_NORM_RELATIVE )))){
      if(!obj_ptr->evaluate(pd, currentOFValue, err)){
        err.set_msg("Invalid patch passed to TerminationCriterion.");
      }
        //std::cout<<"\nOF val "<<currentOFValue;
      MSQ_CHKERR(err);
    }
      //std::cout<<"\nOF current "<<currentOFValue;
      //std::cout<<"\nOF previous - current "<<previousOFValue-currentOFValue;
      // if termination on quality improvement absolute
    if(terminationCriterionFlag & QUALITY_IMPROVEMENT_ABSOLUTE){
        //if the improvement was enough
      if(currentOFValue <= qualityImprovementAbsoluteEps){
        return_flag = true;
      }
    }
     // if termination on quality improvement relative
    if(terminationCriterionFlag & QUALITY_IMPROVEMENT_RELATIVE){
        //if the improvement was enough
      if((currentOFValue-lowerOFBound)<=(qualityImprovementRelativeEps*(
         initialOFValue-lowerOFBound))){
        return_flag = true;
      }
    }
      //if termination on successive improvements absolute
    if(terminationCriterionFlag & SUCCESSIVE_IMPROVEMENTS_ABSOLUTE){
        //if the last improvement was not significant enough
        //PRINT_INFO("\ncur = %f, prev = %f, diff = %f, eps = %f",currentOFValue,previousOFValue, previousOFValue-currentOFValue,successiveImprovementsAbsoluteEps);
      if((previousOFValue-currentOFValue)<=successiveImprovementsAbsoluteEps){
        return_flag = true;
      }
    }
    //if termination on successive improvements relative
    if(terminationCriterionFlag & SUCCESSIVE_IMPROVEMENTS_RELATIVE){
        //if the last improvement was not significant enough
      if((previousOFValue-currentOFValue)<=(successiveImprovementsRelativeEps*(
         initialOFValue-currentOFValue))){
        return_flag = true;
      }
    }
      //if not termination, update the previousOFValue to be currentOFValue
    previousOFValue=currentOFValue;
  }//end of termination criteria which need objective calculated
  
    
    //if terminating on bounded vertex movement (a bounding box for the mesh)
  if(terminationCriterionFlag & BOUNDED_VERTEX_MOVEMENT){
    MsqVertex* vert = pd.get_vertex_array(err);
    int num_vert = pd.num_vertices();
    int i=0;
      //for each vertex
    for(i=0;i<num_vert;++i){
        //if any of the coordinates are greater than eps
      if( (vert[i][0]>boundedVertexMovementEps) ||
          (vert[i][1]>boundedVertexMovementEps) ||
          (vert[i][2]>boundedVertexMovementEps) ){
          //then return true
        return_flag = true;
      }
        //otherwise consider the next vertex
    }
  }
    //if none of the criteria were satisfied
  return return_flag;
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
                                 (GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_INF_NORM_ABSOLUTE)  |
                                 (GRADIENT_L2_NORM_RELATIVE | GRADIENT_INF_NORM_RELATIVE)  |
                                 BOUNDED_VERTEX_MOVEMENT)){
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
#define  __FUNC__ "TerminationCriterion::terminate_with_function_and_gradient"
/*!Sets the function and gradient values to be used in terminate (if needed),
  and then calls terminate using the given PatchData and ObjectiveFunction
  pointer.  Finally resets the functionSupplied and gradientSupplied
  booleans to false for the next iteration.*/
bool TerminationCriterion::terminate_with_function_and_gradient(PatchData &pd, ObjectiveFunction* obj_ptr, double func_val, Vector3D* sup_grad, MsqError &err)
{
  // outputs OF value.
  MSQ_DEBUG_ACTION(1,{std::cout << "  o TermCrit -- OF value: "
                                << func_val << endl;});

  //set functionSupplied and gradientSupplied booleans to true
  functionSupplied=true;
  gradientSupplied=true;
    //set the function value and gradient array
  currentOFValue=func_val;
  suppliedGradientArray=sup_grad;
    //call terminate
  bool return_bool=terminate(pd, obj_ptr,err);
    //reset the booleans to false
  functionSupplied=false;
  gradientSupplied=false;
    //return the value given by terminate
  return return_bool;
}



#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::cull_vertices"
/*!This function checks the culling method criterion supplied to the object
  by the user.  If the user does not supply a culling method criterion,
  the default criterion is NONE, and in that case, no culling is performed.
  If the culling method criterion is satisfied, the interior vertices
  of the given patch are flagged as soft_fixed.  Otherwise, the soft_fixed
  flag is removed from each of the vertices in the patch (interior and
  boundary vertices).  Also, if the criterion was satisfied, then the
  function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::cull_vertices(PatchData &pd,
                                      ObjectiveFunction* obj_ptr,
                                      MsqError &err)
{
    //PRINT_INFO("CULLING_METHOD FLAG = %i",cullingMethodFlag);
  
    //cull_bool will be changed to true if the criterion is satisfied
  bool cull_bool=false;
  switch(cullingMethodFlag){
      //if no culling is requested, always return false
    case NONE:
       return cull_bool;
         //if culling on quality improvement absolute
    case QUALITY_IMPROVEMENT_ABSOLUTE:
         //get objective function value
       if(!obj_ptr->evaluate(pd, currentOFValue, err)){
         err.set_msg("Invalid patch passed to TerminationCriterion.");
         MSQ_CHKERR(err);
       }
         //if the improvement was enough, cull
       if(currentOFValue <= cullingEps)
       {
         cull_bool=true;  
       }
         //PRINT_INFO("\ncurrentOFValue = %f, bool = %i\n",currentOFValue,cull_bool);
       
       break;
         //if culing on quality improvement relative
    case QUALITY_IMPROVEMENT_RELATIVE:
         //get objective function value
       if(!obj_ptr->evaluate(pd, currentOFValue, err)){
         err.set_msg("Invalid patch passed to TerminationCriterion.");
         MSQ_CHKERR(err);
       }
         //if the improvement was enough, cull
       if((currentOFValue-lowerOFBound)<=
          (cullingEps*(initialOFValue-lowerOFBound)))
       {
         cull_bool=true;  
       }
       break;
         //if culling on vertex movement absolute
    case VERTEX_MOVEMENT_ABSOLUTE:
         //if movement was enough, cull
       if(pd.get_max_vertex_movement_squared(previousVerticesMemento,err)<=
          cullingEps){
         cull_bool=true;  
       }
         //PRINT_INFO("\nVertexMovement=-%f, cull =%i",pd.get_max_vertex_movement_squared(previousVerticesMemento,err),cull_bool);
       
       break;
         //if culling on vertex movement relative
    case VERTEX_MOVEMENT_RELATIVE:
         //if movement was small enough, cull
       if(pd.get_max_vertex_movement_squared(previousVerticesMemento,err)<=
          (cullingEps*
           pd.get_max_vertex_movement_squared(initialVerticesMemento,err))){
         cull_bool=true;  
       }
       break;
    default:
       err.set_msg("Requested culling method not yet implemented.");
  };
    //Now actually have patch data cull vertices
  if(cull_bool)
  {
    pd.set_free_vertices_soft_fixed(err);
  }
  else
  {
    pd.set_all_vertices_soft_free(err); 
  }
  MSQ_CHKERR(err);
  return cull_bool;
}

#undef __FUNC__
#define  __FUNC__ "TerminationCriterion::cleanup"
/*!
  Currently this only deletes the memento of the vertex positions and the
  mGrad vector if neccessary.
  When culling, we remove the soft fixed flags from all of the vertices.
 */
void TerminationCriterion::cleanup(MeshSet &ms, MsqError &err)
{
    //clean up the memento if used
  if(totalFlag & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE )){
    delete previousVerticesMemento;
    previousVerticesMemento=NULL;
    if(terminationCriterionFlag & (VERTEX_MOVEMENT_RELATIVE))
    {
      delete initialVerticesMemento;
      initialVerticesMemento=NULL;
    }
  }
  if(totalFlag & ((GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_INF_NORM_ABSOLUTE)
                  | (GRADIENT_L2_NORM_RELATIVE | GRADIENT_INF_NORM_RELATIVE)  ) ){
    delete [] mGrad;
  }
  if(cullingMethodFlag){
    ms.clear_all_soft_fixed_flags(err);
  }
}
