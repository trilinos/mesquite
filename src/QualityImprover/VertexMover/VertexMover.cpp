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

#include <fstream>
#include <iostream>
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
  std::cout << "o Executing VertexMover::loop_over_mesh()\n";
  
    // creates a PatchData object at the VertexMover level
    // in order to reduce the number of memory allocations
  PatchData* patch_data=0;
  bool next_patch=true;
  // if we have already been provided a Global Patch (from a previous algorithm).
  
  if (get_global_patch(err) != 0) {
    if (get_patch_type() != PatchData::GLOBAL_PATCH) {
      err.set_msg("PatchDataUser::globalPatch should be NULL.");
      MSQ_CHKERR(err);
    }
    patch_data = get_global_patch(err);
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

  std::cout<<"\n";
    //skip booleans set to false if the initial mesh (or patch) satisfies
    //the termination criteria.
  bool inner_skip=false;
  bool outer_skip=outer_crit->reset(ms,objFunc,err);
  if(!outer_skip){
    
    while ( !stop_met ) {
        //Status bar
      std::cout<<".";
      std::cout.flush();
        // Prior to looping over the patches.
        // Probably want to pass the MeshSet.  
      this->initialize_mesh_iteration(*patch_data, err);MSQ_CHKERR(err); 

      
      if (get_global_patch(err)==0) {
        // propagates information from QualityImprover to MeshSet
        //try to get the first patch, if no patches can be created
        //skip optimization and terminate.
        std::cout << "get_global_patch== "<< get_global_patch(err) <<" .\n"; //dbg 
        next_patch =  ms.get_next_patch(*patch_data, this, err);
      }
        
      if(!next_patch){
        stop_met=true;
          //PRINT_INFO("\nTerminating due to no more free nodes\n");
          //call terminate() anyway, even though we must terminate
        outer_crit->terminate(ms,objFunc,err);
      }
        //otherwise at least one patch can be created
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
            //if patch is global, don't try to get the next patch
          if (get_patch_type() == PatchData::GLOBAL_PATCH) {
            next_patch = false; }
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
  
  std::cout<<"\n";
    //call the criteria's cleanup funtions.
  outer_crit->cleanup(ms,err);
  inner_crit->cleanup(ms,err);
    //call the optimization cleanup function.
  this->cleanup();
  
  if (get_global_patch(err)==0) {
    delete patch_data;
  }

  return 0.;
}
  
