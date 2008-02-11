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
  \file   SteepestDescent.cpp
  \brief  

  Implements the SteepestDescent class member functions.
  
  \author Thomas Leurent
  \date   2002-06-13
*/

#include "SteepestDescent.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <memory.h>
#else
#  include <memory>
#endif

namespace Mesquite {

msq_std::string SteepestDescent::get_name() const
  { return "SteepestDescent"; }
  
PatchSet* SteepestDescent::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

SteepestDescent::SteepestDescent(ObjectiveFunction* of, bool Nash) 
  : VertexMover(of, Nash),
    PatchSetUser(true)
{
}  
  

void SteepestDescent::initialize(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void SteepestDescent::initialize_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void SteepestDescent::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  MSQ_FUNCTION_TIMER( "SteepestDescent::optimize_vertex_positions" );
    //PRINT_INFO("\no  Performing Steepest Descent optimization.\n");
  // Get the array of vertices of the patch. Free vertices are first.
  const double STEP_DECREASE_FACTOR = 0.5, STEP_INCREASE_FACTOR = 1.1;
  int num_vertices = pd.num_free_vertices();
  msq_std::vector<Vector3D> gradient(num_vertices), dk(num_vertices);
  //double norm;
  bool sd_bool=true;//bool for OF values
  double min_edge_len, max_edge_len;
  double step_size, original_value;
  PatchDataVerticesMemento* pd_previous_coords;
  TerminationCriterion* term_crit=get_inner_termination_criterion();
  OFEvaluator& obj_func = get_objective_function_evaluator();
  
  pd_previous_coords = pd.create_vertices_memento( err ); MSQ_ERRRTN(err);
  msq_std::auto_ptr<PatchDataVerticesMemento> memento_deleter( pd_previous_coords );

  pd.get_minmax_edge_length( min_edge_len, max_edge_len );
  step_size = 10*max_edge_len;

  sd_bool = obj_func.update( pd, original_value, gradient, err ); MSQ_ERRRTN(err);
    //set an error if initial patch is invalid.
  if(!sd_bool){
    MSQ_SETERR(err)("SteepestDescent passed invalid initial patch.",
                    MsqError::INVALID_ARG);
    return;
  }

  // does the steepest descent iteration until stopping is required.
  while (!term_crit->terminate()) {
    MSQ_DBGOUT(3) << "Iteration " << term_crit->get_iteration_count() << msq_stdio::endl;
    MSQ_DBGOUT(3) << "  o  original_value: " << original_value << msq_stdio::endl;

      // save vertex coords
    pd.recreate_vertices_memento( pd_previous_coords, err ); MSQ_ERRRTN(err);

      // computes the gradient norm
    //norm = length( &gradient[0], gradient.size() );
    //MSQ_DBGOUT(3) << "  o  gradient norm: " << norm << msq_stdio::endl;

    if (length_squared(&gradient[0], gradient.size())*step_size*step_size < DBL_EPSILON)
      break;

    // ******** Chooses the search direction ********
    // i.e., -gradient for the steepest descent
    for (int i=0; i<num_vertices; ++i)
      dk[i] = -gradient[i];
    //  for (int j=0; j<3; ++j)
    //    dk[i][j] = -gradient[i][j] / norm;

    // ******* Improve Quality *******
    
      // Loop to find a step size that improves quality
      // (or until max iterations).
    int nb_iter;
    for (nb_iter = 0; nb_iter < 10; ++nb_iter) {
        // change vertices coordinates in PatchData according to descent
        //direction.
      pd.move_free_vertices_constrained(&dk[0], dk.size(), step_size, err);  MSQ_ERRRTN(err);
      // and evaluate the objective function with the new node positions.
      double new_value;
      sd_bool=obj_func.evaluate(pd, new_value, err);  MSQ_ERRRTN(err);

      MSQ_DBGOUT(3) << "    o  step_size: " << step_size << msq_stdio::endl;
      MSQ_DBGOUT(3) << "    o  new_value: " << new_value << msq_stdio::endl;
    
        // If value was reduced by step, then keep it (exit loop)
      if (sd_bool && new_value < original_value) 
        break;
        
      // If here, then try again with smaller step.
      
      // undoes node movement
      pd.set_to_vertices_memento( pd_previous_coords, err );  MSQ_ERRRTN(err);
      // and reduces step size to try again.
      step_size *= STEP_DECREASE_FACTOR;
    }
    if (!sd_bool)
      break;
    
    if (0 == nb_iter)
      step_size *= STEP_INCREASE_FACTOR;
    
    obj_func.update(pd, original_value, gradient, err ); MSQ_ERRRTN(err);
    term_crit->accumulate_inner( pd, original_value, &gradient[0], err ); MSQ_ERRRTN(err); 
    term_crit->accumulate_patch( pd, err );  MSQ_ERRRTN(err);
  }
}


void SteepestDescent::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  //  cout << "- Executing SteepestDescent::iteration_complete()\n";
}
  
void SteepestDescent::cleanup()
{
  //  cout << "- Executing SteepestDescent::iteration_end()\n";
}
  
} // namespace Mesquite
