/*!
  \file   VertexMover.cpp
  \brief  

  The VertexMover Class is the base class for all the smoothing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/


#include "VertexMover.hpp"
#include "MeshSet.hpp"
//#include "StoppingCriterion.hpp"

#include "TSTT_C.h"

using namespace Mesquite;

// VertexMover::VertexMover()
// {
// }


/*! \fn VertexMover::loop_over_mesh

    \brief Improves the quality of the MeshSet, calling some
    methods specified in a class derived from VertexMover

    \param const MeshSet &: this MeshSet is looped over. Only the
    mutable data members are changed (such as currentVertexInd).
  */
#undef __FUNC__
#define __FUNC__ "VertexMover::loop_over_mesh" 
void VertexMover::loop_over_mesh(MeshSet &ms, MsqError &err)
{
  cout << "o Executing VertexMover::loop_over_mesh()\n";
  set_mesh_set(&ms);
  StoppingCriterion* crit = get_stopping_criterion();
  crit->reset_all(err);
  bool stop_met=crit->stop(ms,err); MSQ_CHKERR(err);

  int qi_depth=get_patch_depth();
  int dummy=0;
  ms.set_patch_type(MeshSet::ELEMENTS_ON_VERTEX_PATCH, qi_depth, dummy);
  ms.copy_culling_method_bits( QualityImprover::get_culling_method_bits() );
  
  // creates a PatchData object at the VertexMover level
    // in order to reduce the number of memory allocations
  PatchData patch_data;
  this->initialize(patch_data, err); MSQ_CHKERR(err);
    //reset the stopping criterion's loop counter
  crit->reset_counter();  
  while ( !stop_met ) {
    this->initialize_mesh_iteration(patch_data, err);MSQ_CHKERR(err); 
    while( ms.get_next_patch(patch_data, err) ) {
      MSQ_CHKERR(err);
        this->optimize_vertex_positions(patch_data, err); MSQ_CHKERR(err);
        ms.update_mesh(patch_data, err); MSQ_CHKERR(err); // TSTT mesh update !!
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
  
