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
//#include "StoppingCriterion.hpp"

using namespace Mesquite;

VertexMover::VertexMover() :
  QualityImprover()
{  
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
  set_mesh_set(&ms);
  StoppingCriterion* crit = get_stopping_criterion();
  if(crit==0){
    err.set_msg("Stopping Criterion pointer is Null");
    return;
  }
  crit->reset_all(err);
  bool stop_met=crit->stop(ms,err); MSQ_CHKERR(err);

  // creates a PatchData object at the VertexMover level
    // in order to reduce the number of memory allocations
  PatchData patch_data;
  bool next_patch;
  this->initialize(patch_data, err); MSQ_CHKERR(err);
    //reset the stopping criterion's loop counter
  crit->reset_counter();  
  while ( !stop_met ) {
    this->initialize_mesh_iteration(patch_data, err);MSQ_CHKERR(err); 

    // propagates information from QualityImprover to MeshSet
    int dummy=0;
    int qi_depth=get_patch_depth();
    MeshSet::PatchType patch_type = this->get_patch_type();
    ms.set_patch_type(patch_type, qi_depth, dummy);
    ms.copy_culling_method_bits( QualityImprover::get_culling_method_bits() );
  
    next_patch = true;
    while( next_patch )
    {
      Timer loop_timer;
#if MSQ_DEBUG_LEVEL >= 3
      double aomd_t=0;
      double msq_t=0;
#endif
      MSQ_DEBUG_ACTION(3,{std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
                       loop_timer.since_last_check(); });
      next_patch =  ms.get_next_patch(patch_data, err); 
      MSQ_CHKERR(err);
     
      MSQ_DEBUG_ACTION(3,{aomd_t += loop_timer.since_last_check();
                       std::cout << "\t\t- total time spent retrieving AOMD info: "
                                 << aomd_t << std::endl; });
      if (next_patch == true ) {
        MSQ_DEBUG_ACTION(3,{loop_timer.since_last_check(); });
        this->optimize_vertex_positions(patch_data, err); MSQ_CHKERR(err);
        patch_data.update_mesh(err); MSQ_CHKERR(err); // TSTT mesh update !!
        MSQ_DEBUG_ACTION(3,{msq_t += loop_timer.since_last_check();
                         std::cout << "\t\t- total time optimizing patch: "
                                   << msq_t << std::endl;});
      }
      if (ms.get_patch_type() == MeshSet::GLOBAL_PATCH) {
        next_patch = false; }
    }
    this->terminate_mesh_iteration(patch_data, err); MSQ_CHKERR(err);
      //increment stopping criterion's loop counter
    crit->increment_counter();
      //if global don't loop again.
      //if local loop until criteria are met
    if(qi_depth<1)
      stop_met=0;
    else
      stop_met=crit->stop(ms,err); MSQ_CHKERR(err);
  }
  
  this->cleanup();
}
  
