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

#include <fstream.h>
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
void VertexMover::loop_over_mesh(MeshSet &ms, MsqError &err)
{
  std::cout << "o Executing VertexMover::loop_over_mesh()\n";

  // creates a PatchData object at the VertexMover level
    // in order to reduce the number of memory allocations
  PatchData patch_data;
  bool next_patch;
  TerminationCriterion* outer_crit=this->get_outer_termination_criterion();
  TerminationCriterion* inner_crit=this->get_inner_termination_criterion();
    //if outer-criterion is NULL, we wet an error.
  if(outer_crit == 0){
    err.set_msg("Termination Criterion pointer is Null");
    return;
  }
  outer_crit->initialize(ms, patch_data, err);
    //inner criterion is allowed to be NULL.
  if(inner_crit!=NULL)
    inner_crit->initialize(ms, patch_data, err);
  
  bool stop_met=false; MSQ_CHKERR(err);
  // This should probably pass the MeshSet so that the data requirements
  // can be calculated exactly.
  this->initialize(patch_data, err); MSQ_CHKERR(err);
    //reset the stopping criterion's loop counter
  std::cout<<"\n";
  outer_crit->reset(ms,objFunc,err);
  while ( !stop_met ) {
      //Status bar
    std::cout<<".";
    std::cout.flush();
    // Prior to looping over the patches.
    // Probably want to pass the MeshSet.  
    this->initialize_mesh_iteration(patch_data, err);MSQ_CHKERR(err); 

    // propagates information from QualityImprover to MeshSet
    next_patch = true;
    while( next_patch )
    {
      Timer loop_timer;
#if MSQ_DEBUG_LEVEL != 0
      double aomd_t=0;
      double msq_t=0;
#endif
      MSQ_DEBUG_ACTION(3,{std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
                       loop_timer.since_last_check(); });
      next_patch =  ms.get_next_patch(patch_data, this, err); 
      MSQ_CHKERR(err);
     
      MSQ_DEBUG_ACTION(3,{aomd_t += loop_timer.since_last_check();
                       std::cout << "\t\t- total time spent retrieving AOMD info: "
                                 << aomd_t << std::endl; });
      if (next_patch == true ) {
        MSQ_DEBUG_ACTION(3,{loop_timer.since_last_check(); });
        if(inner_crit!=NULL){
          inner_crit->reset(patch_data,objFunc,err);
        }
        this->optimize_vertex_positions(patch_data, err); MSQ_CHKERR(err);
        if(inner_crit!=NULL){
          inner_crit->cull_vertices(patch_data,objFunc,err);
        }
          //we need update to also update the soft_fixed flags
        patch_data.update_mesh(err); MSQ_CHKERR(err); // TSTT mesh update !!
        MSQ_DEBUG_ACTION(3,{msq_t += loop_timer.since_last_check();
                         std::cout << "\t\t- total time optimizing patch: "
                                   << msq_t << std::endl;});
      }
      if (get_patch_type() == PatchData::GLOBAL_PATCH) {
        next_patch = false; }
    }
    this->terminate_mesh_iteration(patch_data, err); MSQ_CHKERR(err);
      //if global don't loop again.
      //if local loop until criteria are met
    if(this->get_patch_type() == PatchData::GLOBAL_PATCH)
      stop_met=true;
    else
      stop_met=outer_crit->terminate(ms,objFunc,err);
    MSQ_CHKERR(err);
  }
  std::cout<<"\n";
  outer_crit->cleanup(ms,err);
  if(inner_crit!=NULL){ 
    inner_crit->cleanup(ms,err);
  }
  this->cleanup();
}
  
