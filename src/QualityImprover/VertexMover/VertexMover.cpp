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
/*!
  \file   VertexMover.cpp
  \brief  

  The VertexMover Class is the base class for all the smoothing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/


#include "VertexMover.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

namespace Mesquite {

VertexMover::VertexMover() :
  QualityImprover()
{
  objFunc=NULL;
}


/*! \fn VertexMover::loop_over_mesh(MeshSet &ms, MsqError &err)

    \brief Improves the quality of the MeshSet, calling some
    methods specified in a class derived from VertexMover

    \param const MeshSet &: this MeshSet is looped over. Only the
    mutable data members are changed (such as currentVertexInd).
  */
double VertexMover::loop_over_mesh(MeshSet &ms, MsqError &err)
{
  set_mesh_set(&ms);
  
    // creates a PatchData object at the VertexMover level
    // in order to reduce the number of memory allocations
  PatchData local_patch_data;
  PatchData* patch_data=0;
  bool next_patch=true;
  // if we have already been provided a Global Patch (from a previous algorithm).
  
  if (get_global_patch() != 0) {
    if (get_patch_type() != PatchData::GLOBAL_PATCH) {
      MSQ_SETERR(err)("PatchDataUser::globalPatch should be NULL.", MsqError::INVALID_STATE);
      return 0;
    }
    patch_data = get_global_patch();
    next_patch = true; // same as MeshSet::get_next_patch()
  }
  else {
    patch_data = &local_patch_data;
  }
  
  TerminationCriterion* outer_crit=this->get_outer_termination_criterion();
  TerminationCriterion* inner_crit=this->get_inner_termination_criterion();
    //if outer-criterion is NULL, we wet an error.
  if(outer_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer is Null", MsqError::INVALID_STATE);
    return 0.;
  }
    //if inner-criterion is NULL, we wet an error.
  if(inner_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer for inner loop is Null", MsqError::INVALID_STATE);
    return 0.;
  }
  
  bool stop_met = false, inner_skip = false, outer_skip = false;
  
    //initialize both criterion objects
  outer_crit->initialize(ms, *patch_data, err);  if (MSQ_CHKERR(err)) goto ERROR;
  inner_crit->initialize(ms, *patch_data, err);  if (MSQ_CHKERR(err)) goto ERROR;
  
  // This should probably pass the MeshSet so that the data requirements
  // can be calculated exactly.
  this->initialize(*patch_data, err); if (MSQ_CHKERR(err)) goto ERROR;

    //skip booleans set to false if the initial mesh (or patch) satisfies
    //the termination criteria.
  outer_skip=outer_crit->reset(ms,objFunc,err);  
  if (MSQ_CHKERR(err)) goto ERROR;
  if(!outer_skip){
    
    while ( !stop_met ) {
        //Status bar
//      Message::print_info(".");
        // Prior to looping over the patches.
        // Probably want to pass the MeshSet.  
      this->initialize_mesh_iteration(*patch_data, err);
      if (MSQ_CHKERR(err)) goto ERROR;

        // if there is no global patch previously available
      if (get_global_patch()==0) {
        // propagates information from QualityImprover to MeshSet
        //try to get the first patch, if no patches can be created
        //skip optimization and terminate.
        next_patch =  ms.get_next_patch(*patch_data, this, err);
        if (MSQ_CHKERR(err)) goto ERROR;
      }
        
      if(!next_patch){
        stop_met=true;
          //PRINT_INFO("\nTerminating due to no more free nodes\n");
          //call terminate() anyway, even though we must terminate
        outer_crit->terminate(ms,objFunc,err); 
        if (MSQ_CHKERR(err)) goto ERROR;
      }
        //otherwise one patch has been created and more could be created later.
      else{

          //loop over these patches
        while( next_patch )
        {
          if (next_patch == true ) {

            inner_skip=inner_crit->reset(*patch_data,objFunc,err);
              //if inner criteria are initially satisfied, skip opt.
              //otherwise:
            if(!inner_skip){
              this->optimize_vertex_positions(*patch_data, err);  
              if (MSQ_CHKERR(err)) goto ERROR;
            }
            inner_crit->cull_vertices(*patch_data,objFunc,err);
            patch_data->update_mesh(err); // TSTT mesh update !!
            MSQ_ERRZERO(err); 
          }
            // if patch is global, don't try to get the next patch
          if (get_patch_type() == PatchData::GLOBAL_PATCH) {
            next_patch = false; }
            // if patch is local, try to get the next one 
          else{
            next_patch =  ms.get_next_patch(*patch_data, this, err);  
            if (MSQ_CHKERR(err)) goto ERROR;
          }
        }
        this->terminate_mesh_iteration(*patch_data, err); 
        if (MSQ_CHKERR(err)) goto ERROR;
          //check the criteria on the outer loop
        stop_met=outer_crit->terminate(ms,objFunc,err); 
        if (MSQ_CHKERR(err)) goto ERROR;
      }
    } 
  }

ERROR:  
    //call the criteria's cleanup funtions.
  outer_crit->cleanup(ms,err);
  inner_crit->cleanup(ms,err);
    //call the optimization cleanup function.
  this->cleanup();

  return 0.;
}
  
} // namespace Mesquite
