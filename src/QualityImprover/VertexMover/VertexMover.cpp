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
#include "MsqMessage.hpp"

using namespace Mesquite;

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
#undef __FUNC__
#define __FUNC__ "VertexMover::loop_over_mesh" 
double VertexMover::loop_over_mesh(MeshSet &ms, MsqError &err)
{
  set_mesh_set(&ms);
  
    // creates a PatchData object at the VertexMover level
    // in order to reduce the number of memory allocations
  PatchData* patch_data=0;
  bool next_patch=true;
  // if we have already been provided a Global Patch (from a previous algorithm).
  
  if (get_global_patch() != 0) {
    if (get_patch_type() != PatchData::GLOBAL_PATCH) {
      err.set_msg("PatchDataUser::globalPatch should be NULL.");
    }
    patch_data = get_global_patch();
    next_patch = true; // same as MeshSet::get_next_patch()
  }
  else {
    patch_data = new PatchData;
  }
  
  TerminationCriterion* outer_crit=this->get_outer_termination_criterion();
  TerminationCriterion* inner_crit=this->get_inner_termination_criterion();
    //if outer-criterion is NULL, we wet an error.
  if(outer_crit == 0){
    err.set_msg("Termination Criterion pointer is Null");
    return 0.;
  }
    //if inner-criterion is NULL, we wet an error.
  if(inner_crit == 0){
    err.set_msg("Termination Criterion pointer for inner loop is Null");
    return 0.;
  }
    //initialize both criterion objects
  outer_crit->initialize(ms, *patch_data, err);
  inner_crit->initialize(ms, *patch_data, err);
  
  bool stop_met=false; MSQ_CHKERR(err);
  // This should probably pass the MeshSet so that the data requirements
  // can be calculated exactly.
  this->initialize(*patch_data, err); MSQ_CHKERR(err);

  Message::print_info("\n");
    //skip booleans set to false if the initial mesh (or patch) satisfies
    //the termination criteria.
  bool inner_skip=false;
  bool outer_skip=outer_crit->reset(ms,objFunc,err);
  if(!outer_skip){
    
    while ( !stop_met ) {
        //Status bar
      Message::print_info(".");
        // Prior to looping over the patches.
        // Probably want to pass the MeshSet.  
      this->initialize_mesh_iteration(*patch_data, err);MSQ_CHKERR(err); 

        // if there is no global patch previously available
      if (get_global_patch()==0) {
        // propagates information from QualityImprover to MeshSet
        //try to get the first patch, if no patches can be created
        //skip optimization and terminate.
        next_patch =  ms.get_next_patch(*patch_data, this, err);
        MSQ_CHKERR(err);
      }
        
      if(!next_patch){
        stop_met=true;
          //PRINT_INFO("\nTerminating due to no more free nodes\n");
          //call terminate() anyway, even though we must terminate
        outer_crit->terminate(ms,objFunc,err);
      }
        //otherwise one patch has been created and more could be created later.
      else{

          //loop over these patches
        while( next_patch )
        {
          MSQ_CHKERR(err); 
          if (next_patch == true ) {

            inner_skip=inner_crit->reset(*patch_data,objFunc,err);
              //if inner criteria are initially satisfied, skip opt.
              //otherwise:
            if(!inner_skip){
              this->optimize_vertex_positions(*patch_data, err);
              MSQ_CHKERR(err);
            }
            inner_crit->cull_vertices(*patch_data,objFunc,err);
            patch_data->update_mesh(err);// TSTT mesh update !!
            MSQ_CHKERR(err); 
          }
            // if patch is global, don't try to get the next patch
          if (get_patch_type() == PatchData::GLOBAL_PATCH) {
            next_patch = false; }
            // if patch is local, try to get the next one 
          else{
            next_patch =  ms.get_next_patch(*patch_data, this, err);
          }
        }
        this->terminate_mesh_iteration(*patch_data, err); MSQ_CHKERR(err);
          //check the criteria on the outer loop
        stop_met=outer_crit->terminate(ms,objFunc,err);
        MSQ_CHKERR(err);
      }
    } 
  }
  
  Message::print_info("\n");
    //call the criteria's cleanup funtions.
  outer_crit->cleanup(ms,err);
  inner_crit->cleanup(ms,err);
    //call the optimization cleanup function.
  this->cleanup();
  
  if (get_global_patch()==0) {
    delete patch_data;
  }

  return 0.;
}
  
