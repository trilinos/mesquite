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
#include "ParallelHelperInterface.hpp"
#include <algorithm>

namespace MESQUITE_NS {

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
                                    const Settings* settings,
                                    MsqError& err )
{
    // Get the patch data to use for the first iteration
  OFEvaluator& obj_func = get_objective_function_evaluator();
  
  PatchData patch;
  patch.set_mesh( mesh );
  patch.set_domain( domain );
  patch.attach_settings( settings );
  bool one_patch = false, did_some, all_culled;
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
  
  obj_func.initialize( mesh, domain, settings, patch_set, err ); 
  if (MSQ_CHKERR(err)) goto ERROR;
  
  outer_crit->reset_outer(mesh, domain, obj_func, settings, err); 
  if (MSQ_CHKERR(err)) goto ERROR;
  
 
    // if only one patch, get the patch now
  if (one_patch) {
    patch_set->get_patch( patch_list[0], patch_elements, patch_vertices, err );
    if (MSQ_CHKERR(err)) goto ERROR;
    patch.set_mesh_entities( patch_elements, patch_vertices, err );
    if (MSQ_CHKERR(err)) goto ERROR;
  }
  
   // Loop until outer termination criterion is met
  did_some = true;
  while (!outer_crit->terminate())
  {
    if (!did_some) {
      MSQ_SETERR(err)("Inner termiation criterion satisfied for all patches "
                      "without meeting outer termination criterion.  This is "
                      "an infinite loop.  Aborting.", MsqError::INVALID_STATE);
      break;
    }
    did_some = false;
    all_culled = true;

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
        
        if (patch_elements.empty()) { // no more non-culled vertices
          break;
        }
      
        all_culled = false;
        patch.set_mesh_entities( patch_elements, patch_vertices, err );
        if (MSQ_CHKERR(err)) goto ERROR;
      } else {
        ++p_iter;
        all_culled = false;
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
        did_some = true;
        
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
    
    outer_crit->accumulate_outer( mesh, domain, obj_func, settings, err );
    if (MSQ_CHKERR(err)) goto ERROR;
    
    if (all_culled)
      break;
  }


ERROR:  
    //call the criteria's cleanup funtions.
  outer_crit->cleanup(mesh,domain,err);
  inner_crit->cleanup(mesh,domain,err);
    //call the optimization cleanup function.
  this->cleanup();

  return 0.;
}
  

/*! \brief Improves the quality of the MeshSet, calling some
    methods specified in a class derived from VertexMover

    \param const MeshSet &: this MeshSet is looped over. Only the
    mutable data members are changed (such as currentVertexInd).
  */
double VertexMover::loop_over_mesh( ParallelMesh* mesh,
                                    MeshDomain* domain,
                                    const Settings* settings,
                                    MsqError& err )
{
  msq_std::vector<size_t> junk;
  Mesh::VertexHandle vertex_handle;

    // Get the patch data to use for the first iteration
  OFEvaluator& obj_func = get_objective_function_evaluator();
  
  PatchData patch;
  patch.set_mesh( (Mesh*) mesh );
  patch.set_domain( domain );
  patch.attach_settings( settings );

  ParallelHelper* helper = mesh->get_parallel_helper();

  helper->smoothing_init();

  bool did_some, all_culled;
  msq_std::vector<Mesh::VertexHandle> patch_vertices;
  msq_std::vector<Mesh::ElementHandle> patch_elements;
  msq_std::vector<Mesh::VertexHandle> fixed_vertices;
  msq_std::vector<Mesh::VertexHandle> free_vertices;
   
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
  
  PatchSet* patch_set = get_patch_set();
  if (!patch_set) {
    MSQ_SETERR(err)("No PatchSet for QualityImprover!", MsqError::INVALID_STATE);
    return 0.0;
  }
  patch_set->set_mesh( (Mesh*)mesh );

  msq_std::vector<PatchSet::PatchHandle> patch_list;
  patch_set->get_patch_handles( patch_list, err ); MSQ_ERRZERO(err);
  
    // Initialize outer loop
    
  this->initialize(patch, err);        
  if (MSQ_CHKERR(err)) goto ERROR;
  
  obj_func.initialize( (Mesh*)mesh, domain, settings, patch_set, err ); 
  if (MSQ_CHKERR(err)) goto ERROR;
  
  outer_crit->reset_outer( (Mesh*)mesh, domain, obj_func, settings, err); 
  if (MSQ_CHKERR(err)) goto ERROR;
   
   // Loop until outer termination criterion is met
  did_some = true;
  while (!outer_crit->terminate())
  {
    if (!did_some) {
      MSQ_SETERR(err)("Inner termiation criterion satisfied for all patches "
                      "without meeting outer termination criterion.  This is "
                      "an infinite loop.  Aborting.", MsqError::INVALID_STATE);
      break;
    }
    did_some = false;
    all_culled = true;

    ///*** smooth the interior ***////

    // get the fixed vertices (i.e. the ones *not* part of the first independent set)
    helper->compute_first_independent_set(fixed_vertices);
    // sort the fixed vertices
    msq_std::sort(fixed_vertices.begin(), fixed_vertices.end());

      // Loop over each patch
    msq_std::vector<PatchSet::PatchHandle>::iterator p_iter = patch_list.begin();
    while( p_iter != patch_list.end() )
    {
      // loop until we get a non-empty patch.  patch will be empty
      // for culled vertices with element-on-vertex patches
      do {
	patch_set->get_patch( *p_iter, patch_elements, patch_vertices, err );
	if (MSQ_CHKERR(err)) goto ERROR;
	++p_iter;
      } while (patch_elements.empty() && p_iter != patch_list.end()) ;
        
      if (patch_elements.empty()) { // no more non-culled vertices
	break;
      }

      if (patch_vertices.empty()) // global patch hack (means all mesh vertices)
      {
	mesh->get_all_vertices(patch_vertices, err);
      }

      free_vertices.clear();

      for (size_t i = 0; i < patch_vertices.size(); ++i) 
	if (!msq_std::binary_search(fixed_vertices.begin(), fixed_vertices.end(), patch_vertices[i]))
	  free_vertices.push_back(patch_vertices[i]);

      if (free_vertices.empty()) { // all vertices were fixed -> skip patch
	continue;
      }

      all_culled = false;
      patch.set_mesh_entities( patch_elements, free_vertices, err );
      if (MSQ_CHKERR(err)) goto ERROR;
        
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
        did_some = true;
        
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

    helper->communicate_first_independent_set();

    ///*** smooth the boundary ***////

    while (helper->compute_next_independent_set())
    {
      // Loop over all boundary elements
      while(helper->get_next_partition_boundary_vertex(vertex_handle))
      {
	patch_vertices.clear();
	patch_vertices.push_back(vertex_handle);
	patch_elements.clear(); 
	mesh->vertices_get_attached_elements( &vertex_handle, 
                                              1,
                                              patch_elements, 
                                              junk, err );

	all_culled = false;
	patch.set_mesh_entities( patch_elements, patch_vertices, err );

	if (MSQ_CHKERR(err)) goto ERROR;
        
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
	  did_some = true;
	  
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
      helper->communicate_next_independent_set();
    }

    this->terminate_mesh_iteration(patch, err); 
    if (MSQ_CHKERR(err)) goto ERROR;
    
    outer_crit->accumulate_outer( mesh, domain, obj_func, settings, err );
    if (MSQ_CHKERR(err)) goto ERROR;
    
    if (all_culled)
      break;
  }

ERROR: 
    //call the criteria's cleanup funtions.
  outer_crit->cleanup(mesh,domain,err);
  inner_crit->cleanup(mesh,domain,err);
    //call the optimization cleanup function.
  this->cleanup();
    // close the helper
  helper->smoothing_close();

  return 0.;
}

} // namespace Mesquite
