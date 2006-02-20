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
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

namespace Mesquite {

VertexMover::VertexMover() :
  QualityImprover()
{
  objFunc=NULL;
}


/*
  
    +-----------+
    |Reset Outer|
    |Criterion  |
    +-----------+
          |
          V
          +
        /   \
       /Outer\  YES
+--> <Criterion>-----> DONE
|      \Done?/
|       \   /
|         + 
|         |NO
|         V
|   +-----------+
1   |Reset Mesh |
|   | Iteration |
|   +-----------+
|         |
|         V
|         +
|       /  \
|   NO /Next\
+-----<Patch > <-----+
       \    /        |
        \  /         |
          +          |
       YES|          |
          V          |
    +-----------+    |
    |Reset Inner|    |
    |Criterion  |    2
    +-----------+    |
          |          |
          V          |  
          +          |               
        /   \        |
       /Inner\  YES  |
+--> <Criterion>-----+    --------------+
|      \Done?/                          |
|       \   /                           |
|         +                             |
|         |NO                           |
|         V                          Inside
3   +-----------+                    Smoother
|   |   Smooth  |                       |
|   |   Patch   |                       |
|   +-----------+                       |
|         |                             |
----------+               --------------+
                      
*/        
                      


/*! \brief Improves the quality of the MeshSet, calling some
    methods specified in a class derived from VertexMover

    \param const MeshSet &: this MeshSet is looped over. Only the
    mutable data members are changed (such as currentVertexInd).
  */
double VertexMover::loop_over_mesh( Mesh* mesh,
                                    MeshDomain* domain,
                                    PatchData* global_patch,
                                    MsqError& err )
{
    // Get the patch data to use for the first iteration
  bool next_patch = true;
  PatchData local_patch_data;
  PatchData* patch_data = 0;
  const PatchData::PatchType patch_type = get_patch_type();
  
  switch( patch_type )
  {
    default:
      MSQ_SETERR(err)( MsqError::INVALID_ARG,
                       "Invalid patch type for VertexMover (%d).\n", 
                       (int)get_patch_type() );
      return 0;
    
    case PatchData::GLOBAL_PATCH:
      if (global_patch)
      {
        patch_data = global_patch;
      }
      else
      {
        local_patch_data.set_mesh( mesh );
        local_patch_data.set_domain( domain );
        local_patch_data.fill_global_patch( err ); MSQ_ERRZERO(err);
        patch_data = &local_patch_data;
      }
      break;
    
    case PatchData::ELEMENTS_ON_VERTEX_PATCH:
      local_patch_data.set_mesh( mesh );
      local_patch_data.set_domain( domain );
      patch_data = &local_patch_data;
      break;
  }
  
    // Get termination criteria
  TerminationCriterion* outer_crit=this->get_outer_termination_criterion();
  TerminationCriterion* inner_crit=this->get_inner_termination_criterion();
  if(outer_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer is Null", MsqError::INVALID_STATE);
    return 0.;
  }
  if(inner_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer for inner loop is Null", MsqError::INVALID_STATE);
    return 0.;
  }
  
    // If using a local patch, suppress output of inner termination criterion
  if (get_patch_type() != PatchData::GLOBAL_PATCH) {
    inner_crit->set_debug_output_level(3);
  }
  
    // Initialize outer loop
    
  this->initialize(*patch_data, err);        
  if (MSQ_CHKERR(err)) goto ERROR;
  
  outer_crit->reset_outer(mesh, domain, objFunc, err); 
  if (MSQ_CHKERR(err)) goto ERROR;
  
    // Loop until outer termination criterion is met
  while (!outer_crit->terminate())
  {
    if (patch_type == PatchData::ELEMENTS_ON_VERTEX_PATCH)
    {
      size_t free_vtx;
      next_patch = patch_data
        ->get_next_vertex_element_patch( get_nb_layers(), true, free_vtx, err );
      if (MSQ_CHKERR(err)) goto ERROR;
    
      if (!next_patch) // all vertices culled
        break;
    }
  
      // Loop over each patch 
    do
    {
        // Initialize for inner iteration
        
      this->initialize_mesh_iteration(*patch_data, err);
      if (MSQ_CHKERR(err)) goto ERROR;
      
      outer_crit->reset_patch( *patch_data, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
      inner_crit->reset_inner( *patch_data, objFunc, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
      inner_crit->reset_patch( *patch_data, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
        // Don't even call optimizer if inner termination 
        // criterion has already been met.
      if (!inner_crit->terminate())
      {
          // Call optimizer - should loop on inner_crit->terminate()
        this->optimize_vertex_positions( *patch_data, err );
        if (MSQ_CHKERR(err)) goto ERROR;
      
          // Update for changes during inner iteration 
          // (during optimizer loop)
        
        outer_crit->accumulate_patch( *patch_data, err );
        if (MSQ_CHKERR(err)) goto ERROR;
        
        inner_crit->cull_vertices( *patch_data, objFunc, err );
        if (MSQ_CHKERR(err)) goto ERROR;
        
        patch_data->update_mesh( err );
        if (MSQ_CHKERR(err)) goto ERROR;
      }
      
        // If a global patch, then we're done with the inner loop.
      next_patch = false;
      if (patch_type != PatchData::GLOBAL_PATCH)
      {
        size_t free_vtx;
        next_patch = patch_data
          ->get_next_vertex_element_patch( get_nb_layers(), true, free_vtx, err );
        if (MSQ_CHKERR(err)) goto ERROR;
      }
      
    } while (next_patch);

    this->terminate_mesh_iteration(*patch_data, err); 
    if (MSQ_CHKERR(err)) goto ERROR;
    
    if (patch_type != PatchData::GLOBAL_PATCH)
      patch_data->reset_iterators();
    
    outer_crit->accumulate_outer( mesh, domain, err );
    if (MSQ_CHKERR(err)) goto ERROR;
  }


ERROR:  
    //call the criteria's cleanup funtions.
  outer_crit->cleanup(mesh,domain,err);
  inner_crit->cleanup(mesh,domain,err);
    //call the optimization cleanup function.
  this->cleanup();

  return 0.;
}
  
} // namespace Mesquite
