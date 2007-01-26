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
#include "PatchSet.hpp"
#include "PatchData.hpp"

namespace Mesquite {

VertexMover::VertexMover( ObjectiveFunction* OF, bool Nash ) 
  : QualityImprover(),
    objFuncEval( OF, Nash ) 
  {}


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
                                    MappingFunctionSet* map_func,
                                    MsqError& err )
{
    // Get the patch data to use for the first iteration
  OFEvaluator& obj_func = get_objective_function_evaluator();
  
  PatchData patch;
  patch.set_mesh( mesh );
  patch.set_domain( domain );
  patch.set_mapping_functions( map_func );
  bool one_patch = false;
  msq_std::vector<Mesh::VertexHandle> patch_vertices;
  msq_std::vector<Mesh::ElementHandle> patch_elements;
  
  PatchSet* patch_set = get_patch_set();
  if (!patch_set) {
    MSQ_SETERR(err)("No PatchSet for QualityImprover!", MsqError::INVALID_STATE);
    return 0.0;
  }
  patch_set->set_mesh( mesh );
  
  msq_std::vector<PatchSet::PatchHandle> patch_list;
  patch_set->get_patch_handles( patch_list, err ); MSQ_ERRZERO(err);
  
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
  if (patch_list.size() > 1) 
    inner_crit->set_debug_output_level(3);
  else
    one_patch = true;
  
    // Initialize outer loop
    
  this->initialize(patch, err);        
  if (MSQ_CHKERR(err)) goto ERROR;
  
  obj_func.initialize( mesh, domain, map_func, patch_set, err ); 
  if (MSQ_CHKERR(err)) goto ERROR;
  
  outer_crit->reset_outer(mesh, domain, obj_func, err); 
  if (MSQ_CHKERR(err)) goto ERROR;
  
 
    // if only one patch, get the patch now
  if (one_patch) {
    patch_set->get_patch( patch_list[0], patch_elements, patch_vertices, err );
    if (MSQ_CHKERR(err)) goto ERROR;
    patch.set_mesh_entities( patch_elements, patch_vertices, err );
    if (MSQ_CHKERR(err)) goto ERROR;
  }
  
   // Loop until outer termination criterion is met
  while (!outer_crit->terminate())
  {
      // Loop over each patch
    msq_std::vector<PatchSet::PatchHandle>::iterator p_iter = patch_list.begin();
    while( p_iter != patch_list.end() )
    {
      if (!one_patch) { // if only one patch (global) re-use the previous one
          // loop until we get a non-empty patch.  patch will be empty
          // for culled vertices with element-on-vertex patches
        do {
          patch_set->get_patch( *p_iter, patch_elements, patch_vertices, err );
          if (MSQ_CHKERR(err)) goto ERROR;
          ++p_iter;
        } while (patch_elements.empty() && p_iter != patch_list.end()) ;
        
        if (patch_elements.empty()) // no more non-culled vertices
          break;
      
        patch.set_mesh_entities( patch_elements, patch_vertices, err );
        if (MSQ_CHKERR(err)) goto ERROR;
      } else {
        ++p_iter;
      }
        
        // Initialize for inner iteration
        
      this->initialize_mesh_iteration(patch, err);
      if (MSQ_CHKERR(err)) goto ERROR;
      
      obj_func.reset();
      
      outer_crit->reset_patch( patch, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
      inner_crit->reset_inner( patch, obj_func, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
      inner_crit->reset_patch( patch, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
        // Don't even call optimizer if inner termination 
        // criterion has already been met.
      if (!inner_crit->terminate())
      {
          // Call optimizer - should loop on inner_crit->terminate()
        this->optimize_vertex_positions( patch, err );
        if (MSQ_CHKERR(err)) goto ERROR;
      
          // Update for changes during inner iteration 
          // (during optimizer loop)
        
        outer_crit->accumulate_patch( patch, err );
        if (MSQ_CHKERR(err)) goto ERROR;
        
        inner_crit->cull_vertices( patch, obj_func, err );
        if (MSQ_CHKERR(err)) goto ERROR;
        
        patch.update_mesh( err );
        if (MSQ_CHKERR(err)) goto ERROR;
      }
    } 

    this->terminate_mesh_iteration(patch, err); 
    if (MSQ_CHKERR(err)) goto ERROR;
    
    outer_crit->accumulate_outer( mesh, domain, obj_func, err );
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
