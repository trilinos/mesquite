/*!
  \file   Randomize.cpp
  \brief  

  The Randomize Class is the concrete class that randomizes
  the vertex positions.          

  \author Michael Brewer
  \date   2002-10-27
*/


#include "Randomize.hpp"

using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "Randomize::Randomize" 
Randomize::Randomize() 
{
  this->set_name("Randomize");
  mPercent=.05;
}  

#undef __FUNC__
#define __FUNC__ "Randomize::Randomize" 
Randomize::Randomize(double percent) 
{
  this->set_name("Randomize");
  mPercent=percent;
}  
  
#undef __FUNC__
#define __FUNC__ "Randomize::initialize" 
void Randomize::initialize(PatchData &pd, MsqError &err)
{
  this->set_patch_depth(1);
}

#undef __FUNC__
#define __FUNC__ "Randomize::initialize_mesh_iteration" 
void Randomize::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
  //  cout << "- Executing Randomize::iteration_complete()\n";
}

#undef __FUNC__
#define __FUNC__ "Randomize::optimize_vertex_position" 
void Randomize::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
    //cout << "- Executing Randomize::optimize_vertex_position()\n";

  int num_local_vertices = pd.num_vertices();
    //int dim = pd.space_dim();
  int dim = get_mesh_set()->space_dim();  
  // gets the array of coordinates for the patch and print it 
  MsqVertex *patch_coords = pd.get_vertex_array(err); MSQ_CHKERR(err);
  // does the randomize smooth
  randomize_vertex(num_local_vertices-1, &patch_coords[1],
                   patch_coords[0], dim, err); MSQ_CHKERR(err);
}
  
#undef __FUNC__
#define __FUNC__ "Randomize::terminate_mesh_iteration" 
void Randomize::terminate_mesh_iteration(PatchData &pd, MsqError &err)
{
  //  cout << "- Executing Randomize::iteration_complete()\n";
}
  
#undef __FUNC__
#define __FUNC__ "Randomize::cleanup" 
void Randomize::cleanup()
{
  //  cout << "- Executing Randomize::iteration_end()\n";
}
  

