/*!
  \file   TopologyModifier.cpp
  \brief  

  The VertexMover Class is the base class for swapping algorythms etc... 

  \author Thomas Leurent
  \date   2002-01-17
*/


#include "TopologyModifier.hpp"

#include "TSTT_C.h"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "TopologyModifier::loop_over_mesh" 
void TopologyModifier::loop_over_mesh(MeshSet& mesh, MsqError &err)
{
  // TODO: for all vertices
//  cout << "o Executing TopologyModifier::loop_over_mesh()\n"
//       << "  Here we will be looping over elements\n";
  this->iteration_begin();
  this->optimize_connectivity();
  this->iteration_complete();
  this->iteration_end();
}
