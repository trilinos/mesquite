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
#include "MsqVertex.hpp"
#include "MsqInterrupt.hpp"
#include "OFEvaluator.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"
#include "PatchData.hpp"
#include "MeshWriter.hpp"

#include <sstream>

namespace Mesquite {


const unsigned long GRAD_FLAGS = TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE |
                                 TerminationCriterion::GRADIENT_INF_NORM_ABSOLUTE |
                                 TerminationCriterion::GRADIENT_L2_NORM_RELATIVE |
                                 TerminationCriterion::GRADIENT_INF_NORM_RELATIVE;
const unsigned long OF_FLAGS   = TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE |
                                 TerminationCriterion::QUALITY_IMPROVEMENT_RELATIVE |
                                 TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_ABSOLUTE |
                                 TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE;

/*!Constructor initializes all of the data members which are not
  necessarily automatically initialized in their constructors.*/
TerminationCriterion::TerminationCriterion() 
  : mGrad(8),
    initialVerticesMemento(0),
    previousVerticesMemento(0),
    debugLevel(2),
    timeStepFileType(NOTYPE)
{
  terminationCriterionFlag=NONE;
  cullingMethodFlag=NONE;
  cullingEps=0.0;
  initialOFValue=0.0;
  previousOFValue=0.0;
  currentOFValue = 0.0;
  lowerOFBound=0.0;
  initialGradL2NormSquared=0.0;
  initialGradInfNorm=0.0;
    //initial size of the gradient array
  gradL2NormAbsoluteEpsSquared=0.0;
  gradL2NormRelativeEpsSquared=0.0;
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
  
}



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
       gradL2NormAbsoluteEpsSquared=eps*eps;
       break; 
    case GRADIENT_INF_NORM_ABSOLUTE:
       terminationCriterionFlag|=GRADIENT_INF_NORM_ABSOLUTE;
       gradInfNormAbsoluteEps=eps;
       break;
    case GRADIENT_L2_NORM_RELATIVE:
       terminationCriterionFlag|=GRADIENT_L2_NORM_RELATIVE;
       gradL2NormRelativeEpsSquared=eps*eps;
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
       MSQ_SETERR(err)("TCType not valid for this function.",MsqError::INVALID_ARG);
  };
}

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
       MSQ_SETERR(err)("TCType not valid for this function.",MsqError::INVALID_ARG);
  };
}
/*! Function to remove a previously set criterion type.*/
void TerminationCriterion::remove_criterion_type(TCType tc_type,
                                                 MsqError &/*err*/)
{
  terminationCriterionFlag&=(~tc_type);
}


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
       MSQ_SETERR(err)("TCType not valid for this function.",MsqError::INVALID_ARG);
  };
  cullingEps=eps;
}

/*!Sets the culling type to be NONE.*/
void TerminationCriterion::remove_culling(MsqError &/*err*/)
{
  cullingMethodFlag=NONE;
}

  

/*!This version of reset is called using a MeshSet, which implies
  it is only called when this criterion is used as the 'outer' termination
  criterion.  
 */
void TerminationCriterion::reset_outer(Mesh* mesh, 
                                       MeshDomain* domain,
                                       OFEvaluator& obj_eval,
                                       MappingFunctionSet* map_func,
                                       MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  PatchData global_patch;
  global_patch.set_mapping_functions( map_func );
  
    //if we need to fill out the global patch data object.
  if ((totalFlag & (GRAD_FLAGS | OF_FLAGS | VERTEX_MOVEMENT_RELATIVE))
     || timeStepFileType)
  {
    global_patch.set_mesh( mesh );
    global_patch.set_domain( domain );
    global_patch.fill_global_patch( err );
    MSQ_ERRRTN(err);
  }

    //now call the other reset
  reset_inner( global_patch, obj_eval, err ); MSQ_ERRRTN(err);
}
    
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
void TerminationCriterion::reset_inner(PatchData &pd, OFEvaluator& obj_eval,
                                    MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  
    // clear flag for BOUNDED_VERTEX_MOVEMENT
  vertexMovementExceedsBound = 0;
  
    // Use -1 to denote that this isn't initialized yet.
    // As all valid values must be >= 0.0, a negative
    // value indicates that it is uninitialized and is
    // always less than any valid value.
  maxSquaredMovement = -1;
  
    // Clear the iteration count.
  iterationCounter = 0;
  
    //reset the inner timer if needed
  if(totalFlag & CPU_TIME){
    mTimer.reset();
  }
   
    //GRADIENT
  currentGradInfNorm = initialGradInfNorm = 0.0;
  currentGradL2NormSquared = initialGradL2NormSquared = 0.0;
  if(totalFlag & GRAD_FLAGS)
  {
    if (!obj_eval.have_objective_function()) {
      MSQ_SETERR(err)("Error termination criteria set which uses objective "
                      "functions, but no objective function is available.",
                      MsqError::INVALID_STATE);   
      return;
    } 
    int num_vertices=pd.num_free_vertices();
    mGrad.resize( num_vertices );

      //get gradient and make sure it is valid
    bool b = obj_eval.evaluate(pd, currentOFValue, mGrad, err); MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Initial patch is invalid for gradient computation.", 
                      MsqError::INVALID_STATE);
      return;
    } 

      //get the gradient norms
    if (totalFlag & (GRADIENT_INF_NORM_ABSOLUTE|GRADIENT_INF_NORM_RELATIVE))
    {
      currentGradInfNorm = initialGradInfNorm = Linf(&mGrad[0], num_vertices);
      MSQ_DBGOUT(debugLevel) << "  o Initial gradient Inf norm: " 
        << initialGradInfNorm << msq_stdio::endl;
    }  
      
    if (totalFlag & (GRADIENT_L2_NORM_ABSOLUTE|GRADIENT_L2_NORM_RELATIVE))
    {
      currentGradL2NormSquared = initialGradL2NormSquared = length_squared(&mGrad[0], num_vertices);
      MSQ_DBGOUT(debugLevel) << "  o Initial gradient L2 norm: " 
        << msq_stdc::sqrt(initialGradL2NormSquared) << msq_stdio::endl;
    }  

      //the OFvalue comes for free, so save it
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
  }
  //find the initial objective function value if needed and not already
  //computed.  If we needed the gradient, we have the OF value for free.
  // Also, if possible, get initial OF value if writing plot file.  Solvers
  // often supply the OF value for subsequent iterations so by calculating
  // the initial value we can generate OF value plots.
  else if ((totalFlag & OF_FLAGS) || 
           (plotFile.is_open() && pd.num_free_vertices() && obj_eval.have_objective_function()))
  {
      //ensure the obj_ptr is not null
    if(!obj_eval.have_objective_function()){
      MSQ_SETERR(err)("Error termination criteria set which uses objective "
                      "functions, but no objective function is available.",
                      MsqError::INVALID_STATE);
      return;
    }
    
    bool b = obj_eval.evaluate(pd, currentOFValue, err); MSQ_ERRRTN(err);
    if (!b){
      MSQ_SETERR(err)("Initial patch is invalid for evaluation.",MsqError::INVALID_STATE);
      return;
    }
      //std::cout<<"\nReseting initial of value = "<<initialOFValue;
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
  }
  
  if (totalFlag & (GRAD_FLAGS|OF_FLAGS))
    MSQ_DBGOUT(debugLevel) << "  o Initial OF value: " << initialOFValue << msq_stdio::endl;
  
    // Store current vertex locations now, because we'll
    // need them later to compare the current movement with.
  if (totalFlag & VERTEX_MOVEMENT_RELATIVE)
  {
    if (initialVerticesMemento)
    {
      pd.recreate_vertices_memento( initialVerticesMemento, err );
    }
    else
    {
      initialVerticesMemento = pd.create_vertices_memento( err );
    }
    MSQ_ERRRTN(err);
    maxSquaredInitialMovement = DBL_MAX;
  }
  else {
    maxSquaredInitialMovement = 0;
  }

  if (timeStepFileType)
    write_timestep( pd, err);
    
  if (plotFile.is_open()) {
      // two newlines so GNU plot knows that we are starting a new data set
    plotFile << msq_stdio::endl << msq_stdio::endl;
      // write column headings as comment in data file
    plotFile << "#Iter\tCPU\tObjFunc\tGradL2\tGradInf\tMovement" << msq_stdio::endl;
      // write initial values
    plotFile << 0 
     << '\t' << mTimer.since_birth() 
     << '\t' << initialOFValue 
     << '\t' << msq_std::sqrt( currentGradL2NormSquared ) 
     << '\t' << currentGradInfNorm 
     << '\t' << 0.0
     << msq_stdio::endl;
  }
}

void TerminationCriterion::reset_patch(PatchData &pd, MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  if (totalFlag & (VERTEX_MOVEMENT_ABSOLUTE | VERTEX_MOVEMENT_RELATIVE))
  {
    if (previousVerticesMemento)
      pd.recreate_vertices_memento(previousVerticesMemento,err); 
    else
      previousVerticesMemento = pd.create_vertices_memento(err);
    MSQ_ERRRTN(err);
  }
}

void TerminationCriterion::accumulate_inner( PatchData& pd, 
                                             OFEvaluator& of_eval,
                                             MsqError& err )
{
  double of_value = 0;
  
  if (terminationCriterionFlag & GRAD_FLAGS)
  {
    mGrad.resize( pd.num_free_vertices() );
    bool b = of_eval.evaluate(pd, of_value, mGrad, err);
    MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Initial patch is invalid for gradient compuation.",
                      MsqError::INVALID_MESH);
      return;
    }
  }
  else if (terminationCriterionFlag & OF_FLAGS)
  {
    bool b = of_eval.evaluate(pd, of_value, err); MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Invalid patch passed to TerminationCriterion.",
                      MsqError::INVALID_MESH);
      return;
    }
  }

  accumulate_inner( pd, of_value, &mGrad[0], err );  MSQ_CHKERR(err);
}


void TerminationCriterion::accumulate_inner( PatchData& pd, 
                                             double of_value,
                                             Vector3D* grad_array,
                                             MsqError& err )
{
  //if terminating on the norm of the gradient
  //currentGradL2NormSquared = HUGE_VAL;
  if (terminationCriterionFlag & (GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_L2_NORM_RELATIVE)) 
  {
    currentGradL2NormSquared = length_squared(grad_array, pd.num_free_vertices()); // get the L2 norm
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- gradient L2 norm: " 
      << msq_stdc::sqrt(currentGradL2NormSquared) << msq_stdio::endl;
  }
  //currentGradInfNorm = 10e6;
  if (terminationCriterionFlag & (GRADIENT_INF_NORM_ABSOLUTE | GRADIENT_INF_NORM_RELATIVE)) 
  {
    currentGradInfNorm = Linf(grad_array, pd.num_free_vertices()); // get the Linf norm
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- gradient Inf norm: " 
      << currentGradInfNorm << msq_stdio::endl;
  } 
  
  if (terminationCriterionFlag & VERTEX_MOVEMENT_RELATIVE)
  {
    maxSquaredInitialMovement = pd.get_max_vertex_movement_squared(
                               initialVerticesMemento, err );  MSQ_ERRRTN(err);
  }
  
  previousOFValue = currentOFValue;
  currentOFValue = of_value;
  if (terminationCriterionFlag & OF_FLAGS) {
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- OF Value: " << of_value << msq_stdio::endl;
  }
  else if (grad_array) {
    MSQ_DBGOUT(debugLevel) << "  o OF Value: " << of_value << msq_stdio::endl;
  }
  
  ++iterationCounter;
  if (timeStepFileType)
    write_timestep( pd, err);
    
  if (plotFile.is_open()) 
    plotFile << iterationCounter 
     << '\t' << mTimer.since_birth() 
     << '\t' << of_value 
     << '\t' << msq_std::sqrt( currentGradL2NormSquared ) 
     << '\t' << currentGradInfNorm 
     << '\t' << (maxSquaredMovement > 0.0 ? msq_std::sqrt( maxSquaredMovement ) : 0.0)
     << msq_stdio::endl;
}


void TerminationCriterion::accumulate_outer(Mesh* mesh, 
                                            MeshDomain* domain, 
                                            OFEvaluator& of_eval,
                                            MappingFunctionSet* mfs,
                                            MsqError &err)
{
  PatchData global_patch;
  global_patch.set_mapping_functions( mfs );
  
    //if we need to fill out the global patch data object.
  if ((terminationCriterionFlag & (GRAD_FLAGS|OF_FLAGS|VERTEX_MOVEMENT_RELATIVE))
      || timeStepFileType)
  {
    global_patch.set_mesh( mesh );
    global_patch.set_domain( domain );
    global_patch.fill_global_patch( err ); MSQ_ERRRTN(err);
  }
  
  accumulate_inner( global_patch, of_eval, err ); MSQ_ERRRTN(err);
}


void TerminationCriterion::accumulate_patch( PatchData& pd, MsqError& err )
{
  if (terminationCriterionFlag & (VERTEX_MOVEMENT_ABSOLUTE|VERTEX_MOVEMENT_RELATIVE))
  {
    double patch_max_dist = pd.get_max_vertex_movement_squared( previousVerticesMemento, err );
    if (patch_max_dist > maxSquaredMovement)
      maxSquaredMovement = patch_max_dist;
    pd.recreate_vertices_memento( previousVerticesMemento, err );  MSQ_ERRRTN(err);
  }
    
    //if terminating on bounded vertex movement (a bounding box for the mesh)
  if(terminationCriterionFlag & BOUNDED_VERTEX_MOVEMENT)
  {
    MsqVertex* vert = pd.get_vertex_array(err);
    int num_vert = pd.num_free_vertices();
    int i=0;
      //for each vertex
    for(i=0;i<num_vert;++i)
    {
        //if any of the coordinates are greater than eps
      if( (vert[i][0]>boundedVertexMovementEps) ||
          (vert[i][1]>boundedVertexMovementEps) ||
          (vert[i][2]>boundedVertexMovementEps) )
      {
        ++vertexMovementExceedsBound;
      }
    }
  }
}


/*!  This function evaluates the needed information and then evaluates
  the termination criteria.  If any of the selected criteria are satisfied,
  the function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::terminate( )
{
  bool return_flag = false;
  //  cout<<"\nInside terminate(pd,of,err):  flag = "<<terminationCriterionFlag << endl;

    //First check for an interrupt signal
  if (MsqInterrupt::interrupt())
  {
     MSQ_DBGOUT(debugLevel) << "  o TermCrit -- INTERRUPTED" << msq_stdio::endl;
    return true;
  }
  
    //if terminating on numbering of inner iterations
  if (NUMBER_OF_ITERATES & terminationCriterionFlag
    && iterationCounter >= iterationBound)
  {
    return_flag = true;
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- Reached " << iterationBound << " iterations." << msq_stdio::endl;
  }
  
  if (CPU_TIME & terminationCriterionFlag && mTimer.since_birth()>=timeBound)
  {
    return_flag=true;
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- Exceeded CPU time." << msq_stdio::endl;
  }
  
  
  if ((VERTEX_MOVEMENT_ABSOLUTE|VERTEX_MOVEMENT_RELATIVE) & terminationCriterionFlag
      && maxSquaredMovement >= 0.0)
  {
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- Maximuim vertex movement: "
    << sqrt(maxSquaredMovement) << msq_stdio::endl;

    if (VERTEX_MOVEMENT_ABSOLUTE & terminationCriterionFlag 
        && maxSquaredMovement <= vertexMovementAbsoluteEps)
    {
      return_flag = true;
    }
  
    if (VERTEX_MOVEMENT_RELATIVE & terminationCriterionFlag
        && maxSquaredMovement <= vertexMovementRelativeEps*maxSquaredInitialMovement)
    {
      return_flag = true;
    }
  }

  if (GRADIENT_L2_NORM_ABSOLUTE & terminationCriterionFlag &&
      currentGradL2NormSquared <= gradL2NormAbsoluteEpsSquared)
  {
    return_flag = true;
  }
  
  if (GRADIENT_INF_NORM_ABSOLUTE & terminationCriterionFlag &&
      currentGradInfNorm <= gradInfNormAbsoluteEps)
  {
    return_flag = true;
  }
  
  if (GRADIENT_L2_NORM_RELATIVE & terminationCriterionFlag &&
      currentGradL2NormSquared <= (gradL2NormRelativeEpsSquared * initialGradL2NormSquared))
  {
    return_flag = true;
  }
  
  if (GRADIENT_INF_NORM_RELATIVE & terminationCriterionFlag &&
      currentGradInfNorm <= (gradInfNormRelativeEps * initialGradInfNorm))
  {
    return_flag = true;
  }
    //Quality Improvement and Successive Improvements are below.
    // The relative forms are only valid after the first iteration.
  if (QUALITY_IMPROVEMENT_ABSOLUTE & terminationCriterionFlag &&
      currentOFValue <= qualityImprovementAbsoluteEps)
  {
    return_flag = true;
  }
  
    //only valid if an iteration has occurred, see above.
  if(iterationCounter > 0){
    if (SUCCESSIVE_IMPROVEMENTS_ABSOLUTE & terminationCriterionFlag &&
        (previousOFValue - currentOFValue) <= successiveImprovementsAbsoluteEps)
    {  
      return_flag = true;
    }
    if (QUALITY_IMPROVEMENT_RELATIVE & terminationCriterionFlag &&
        (currentOFValue - lowerOFBound) <= 
        qualityImprovementRelativeEps * (initialOFValue - lowerOFBound))
    {
      return_flag = true;
    }
    if (SUCCESSIVE_IMPROVEMENTS_RELATIVE & terminationCriterionFlag &&
        (previousOFValue - currentOFValue) <= 
        successiveImprovementsRelativeEps * (initialOFValue - currentOFValue))
    {
      return_flag = true;
    }
  }
  
  if (BOUNDED_VERTEX_MOVEMENT & terminationCriterionFlag && vertexMovementExceedsBound)
  {
    return_flag = true;
    MSQ_DBGOUT(debugLevel) << "  o TermCrit -- " << vertexMovementExceedsBound
                           << " vertices out of bounds." << msq_stdio::endl;
  }
  
    // clear this value at the end of each iteration
  vertexMovementExceedsBound = 0;
  maxSquaredMovement = -1.0;

  if (timeStepFileType == GNUPLOT && return_flag) {
    MsqError err;
    MeshWriter::write_gnuplot_overlay( iterationCounter, timeStepFileName.c_str(), err );
  }
  
    //if none of the criteria were satisfied
  return return_flag;
}


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
                                      OFEvaluator& of_eval,
                                      MsqError &err)
{
    //PRINT_INFO("CULLING_METHOD FLAG = %i",cullingMethodFlag);
  
    //cull_bool will be changed to true if the criterion is satisfied
  bool b, cull_bool=false;
  double prev_m, init_m;
  switch(cullingMethodFlag){
      //if no culling is requested, always return false
    case NONE:
       return cull_bool;
         //if culling on quality improvement absolute
    case QUALITY_IMPROVEMENT_ABSOLUTE:
         //get objective function value
       b = of_eval.evaluate(pd, currentOFValue, err);
       if (MSQ_CHKERR(err)) return false;
       if (!b) {
         MSQ_SETERR(err)(MsqError::INVALID_MESH);
         return false;
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
       b = of_eval.evaluate(pd, currentOFValue, err);
       if (MSQ_CHKERR(err)) return false;
       if(!b){
         MSQ_SETERR(err)(MsqError::INVALID_MESH);
         return false;
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
       prev_m = pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
       MSQ_ERRZERO(err);
       if(prev_m <= cullingEps){
         cull_bool=true;  
       }
       
       break;
         //if culling on vertex movement relative
    case VERTEX_MOVEMENT_RELATIVE:
         //if movement was small enough, cull
       prev_m = pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
       MSQ_ERRZERO(err);
       init_m = pd.get_max_vertex_movement_squared(initialVerticesMemento,err);
       MSQ_ERRZERO(err);
       if(prev_m <= (cullingEps * init_m)){
         cull_bool=true;  
       }
       break;
    default:
       MSQ_SETERR(err)("Requested culling method not yet implemented.",
                       MsqError::NOT_IMPLEMENTED);
       return false;
  };
    //Now actually have patch data cull vertices
  if(cull_bool)
  {
    pd.set_free_vertices_soft_fixed(err); MSQ_ERRZERO(err);
  }
  else
  {
    pd.set_all_vertices_soft_free(err); MSQ_ERRZERO(err);
  }
  return cull_bool;
}

/*!
  Currently this only deletes the memento of the vertex positions and the
  mGrad vector if neccessary.
  When culling, we remove the soft fixed flags from all of the vertices.
 */
void TerminationCriterion::cleanup(Mesh* mesh, MeshDomain*, MsqError &err)
{
  delete previousVerticesMemento;
  delete initialVerticesMemento;
  previousVerticesMemento = 0;
  initialVerticesMemento = 0;
  
  if (!cullingMethodFlag)
    return;
  
    // Clear soft fixed flag on all vertices
  std::vector<Mesh::VertexHandle> vertices;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRRTN(err);
  std::vector<unsigned char> zeros(vertices.size(), 0);
  mesh->vertices_set_byte( &vertices[0], &zeros[0], vertices.size(), err );
  MSQ_ERRRTN(err);
}

void TerminationCriterion::write_timestep( PatchData& pd, MsqError& err )
{
  pd.update_mesh(err); MSQ_ERRRTN(err);
  std::ostringstream str;
  if (timeStepFileType == VTK) {
    str << timeStepFileName << '_' << iterationCounter << ".vtk";
    MeshWriter::write_vtk( pd.get_mesh(), str.str().c_str(), err );
  }
  else if (timeStepFileType == GNUPLOT) {
    str << timeStepFileName << '.' << iterationCounter;
    MeshWriter::write_gnuplot( pd.get_mesh(), str.str().c_str(), err );
  }
}

void TerminationCriterion::write_iterations( const char* filename, MsqError& err )
{
  if (filename) {
    plotFile.open( filename, msq_stdio::ios::trunc );
    if (!plotFile)
      MSQ_SETERR(err)( MsqError::FILE_ACCESS, "Failed to open plot data file: '%s'", filename );
  }
  else {
    plotFile.close();
  }
}
    


} //namespace Mesquite
