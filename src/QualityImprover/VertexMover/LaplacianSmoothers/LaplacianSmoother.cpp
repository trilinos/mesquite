/*!
  \file   LaplacianSmoother.cpp
  \brief  

  The LaplacianSmoother Class is the concrete class
  that performs Laplacian Smoothing

  \author Thomas Leurent
  \date   2002-01-17
*/

#include "LaplacianSmoother.hpp"

#include "TSTT_C.h"

using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::LaplacianSmoother" 
LaplacianSmoother::LaplacianSmoother() 
{
  this->set_name("LaplacianSmoother");
}  
  
  
#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::initialize" 
void LaplacianSmoother::initialize(PatchData &pd, MsqError &err)
{
  this->set_patch_depth(1);
}

#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::initialize_mesh_iteration" 
void LaplacianSmoother::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
  //  cout << "- Executing LaplacianSmoother::iteration_complete()\n";
}

#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::optimize_vertex_position" 
void LaplacianSmoother::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  std::cout << "- Executing LaplacianSmoother::optimize_vertex_position()\n";

  int num_local_vertices = pd.num_vertices();
    //int dim = pd.space_dim();
  int dim = get_mesh_set()->space_dim();
  
  // gets the array of coordinates for the patch and print it 
  MsqVertex *patch_coords = pd.get_vertex_array(err); MSQ_CHKERR(err);
    //for (size_t i=0; i<num_local_vertices; i++) 
      //cout << "vertex " << i << " : " << patch_coords[i];
  // does the dumb Laplacian smoothing
  centroid_smooth_mesh(num_local_vertices-1, &patch_coords[1],
                       patch_coords[0], dim, err); MSQ_CHKERR(err);
}
  
#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::terminate_mesh_iteration" 
void LaplacianSmoother::terminate_mesh_iteration(PatchData &pd, MsqError &err)
{
  //  cout << "- Executing LaplacianSmoother::iteration_complete()\n";
}
  
#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::cleanup" 
void LaplacianSmoother::cleanup()
{
  //  cout << "- Executing LaplacianSmoother::iteration_end()\n";
}
  

