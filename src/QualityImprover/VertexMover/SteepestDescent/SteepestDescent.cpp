/*!
  \file   SteepestDescent.cpp
  \brief  

  Implements the SteepestDescent class member functions.
  
  \author Thomas Leurent
  \date   2002-06-13
*/

#include "SteepestDescent.hpp"

//Michael delete
#include "AspectRatioGammaQualityMetric.hpp"
#include <time.h>
#include "MsqMessage.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "SteepestDescent::SteepestDescent" 
SteepestDescent::SteepestDescent(ObjectiveFunction* of) :
  VertexMover(),
  objFunc(of)
{
  MsqError err;
  gradientLessThan=.01;
  maxIteration=6;
  this->set_name("SteepestDescent");
  set_patch_type(PatchData::GLOBAL_PATCH, err);
}  
  

#undef __FUNC__
#define __FUNC__ "SteepestDescent::initialize" 
void SteepestDescent::initialize(PatchData &pd, MsqError &err)
{
}

#undef __FUNC__
#define __FUNC__ "SteepestDescent::initialize_mesh_iteration" 
void SteepestDescent::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
}

#undef __FUNC__
#define __FUNC__ "SteepestDescent::optimize_vertex_positions" 
void SteepestDescent::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  PRINT_INFO("\no  Performing Steepest Descent optimization.\n");
  // Get the array of vertices of the patch. Free vertices are first.
  MeshSet *vertex_mover_mesh=get_mesh_set();
  int num_free_vertices = pd.num_free_vertices();
  Vector3D* gradient = new Vector3D[num_free_vertices];
  Vector3D* dk = new Vector3D[num_free_vertices];
  int nb_iterations = 0;
  double norm=10e6;

  double smallest_edge = 0.4; // TODO -- update -- used for step_size
  bool inner_criterion=inner_criterion_met(*vertex_mover_mesh,err);

  objFunc->set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
  
  // does the steepest descent iteration until stopping is required.
  while ( (nb_iterations<maxIteration &&
          norm>gradientLessThan ) && !inner_criterion) {
    
    ++nb_iterations;
  
    objFunc->compute_gradient(pd, gradient, err, num_free_vertices); MSQ_CHKERR(err);

    // Prints out free vertices coordinates. 
    MSQ_DEBUG_ACTION(3,{
      std::cout << "\n  o Free vertices ("<< num_free_vertices <<")original coordinates:\n ";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_CHKERR(err);
      for (int i=0; i<num_free_vertices; ++i)
        std::cout << "\t\t\t" << toto1[i];
    });
      
      // computes the gradient norm
    norm=0;
    for (int n=0; n<num_free_vertices; ++n) 
      norm += gradient[n] % gradient[n]; // dot product
    norm = sqrt(norm);
    MSQ_DEBUG_ACTION(3,{std::cout<< "  o  gradient norm: " << norm << std::endl;});
  
    // ******** Chooses the search direction ********
    // i.e., -gradient for the steepest descent
    for (int i=0; i<num_free_vertices; ++i)
      for (int j=0; j<3; ++j)
        dk[i][j] = -gradient[i][j] / norm;

    // ******* Improve Quality *******
    double original_value = objFunc->evaluate(pd, err);  MSQ_CHKERR(err);
    MSQ_DEBUG_ACTION(3,{std::cout << "  o  original_value: " << original_value << std::endl;});
    
    double new_value = original_value+1;
    // reduces the step size until we get an improvement
    int nb_iter = 0;

    // saves the PatchData coordinates in a memento
    PatchDataCoordsMemento* pd_previous_coords;
    pd_previous_coords = pd.create_coords_memento(err); MSQ_CHKERR(err);
    // Loop to find a step size that improves quality
    double step_size = smallest_edge;
    while (new_value > original_value
           && nb_iter < 10 ) {
      nb_iter++;
      // change vertices coordinates in PatchData according to descent direction.
      pd.move_free_vertices(dk, num_free_vertices, step_size, err); MSQ_CHKERR(err);
      // and evaluate the objective function with the new node positions.
      new_value = objFunc->evaluate(pd, err); MSQ_CHKERR(err);
      MSQ_DEBUG_ACTION(3,{std::cout << "    o  step_size: " << step_size << std::endl; std::cout << "    o  new_value: " << new_value << std::endl;
      });
      
      // if no improvement
      if (new_value > original_value) {
        // undoes node movement
        pd.set_to_coords_memento( pd_previous_coords, err ); MSQ_CHKERR(err);
        // and reduces step size to try again.
        step_size /= 2;
      }
    }

    // Prints out free vertices coordinates. 
    MSQ_DEBUG_ACTION(3,{
      std::cout << "  o Free vertices new coordinates: \n";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_CHKERR(err);
      for (int i=0; i<num_free_vertices; ++i)
        std::cout << "\t\t\t" << toto1[i];
    });
    
    delete pd_previous_coords; // user manages the memento.
    inner_criterion=inner_criterion_met(*vertex_mover_mesh,err); MSQ_CHKERR(err);
    
  }

  delete[] gradient;
  delete[] dk;
    
}


#undef __FUNC__
#define __FUNC__ "SteepestDescent::terminate_mesh_iteration" 
void SteepestDescent::terminate_mesh_iteration(PatchData &pd, MsqError &err)
{
  //  std::cout << "- Executing SteepestDescent::iteration_complete()\n";
}
  
#undef __FUNC__
#define __FUNC__ "SteepestDescent::cleanup" 
void SteepestDescent::cleanup()
{
  //  std::cout << "- Executing SteepestDescent::iteration_end()\n";
}
  

