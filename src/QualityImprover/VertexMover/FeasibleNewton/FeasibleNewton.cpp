// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 15-Jan-03 at 08:05:56
//  LAST-MOD: 16-Jan-03 at 17:04:41 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*!
  \file   FeasibleNewton.cpp
  \brief  

  Implements the FeasibleNewton class member functions.
  
  \author Thomas Leurent
  \date   2003-01-15
*/
// DESCRIP-END.
//

#include "FeasibleNewton.hpp"
#include "MsqFreeVertexIndexIterator.hpp"

using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "FeasibleNewton::FeasibleNewton" 
FeasibleNewton::FeasibleNewton(ObjectiveFunction* of) :
  VertexMover(),
  objFunc(of)
{
  MsqError err;
  gradientLessThan=.01;
  maxIteration=6;
  this->set_name("FeasibleNewton");
  set_patch_type(PatchData::GLOBAL_PATCH, err);
}  
  

#undef __FUNC__
#define __FUNC__ "FeasibleNewton::initialize" 
void FeasibleNewton::initialize(PatchData &pd, MsqError &err)
{
}

#undef __FUNC__
#define __FUNC__ "FeasibleNewton::initialize_mesh_iteration" 
void FeasibleNewton::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
}

#undef __FUNC__
#define __FUNC__ "FeasibleNewton::optimize_vertex_positions" 
void FeasibleNewton::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  PRINT_INFO("\no  Performing Feasible Newton optimization.\n");
  int num_free_vertices = pd.num_free_vertices(err);
  int nb_iterations = 0;
  bool fn_bool=true;// bool used for determining validity of patch
  /* Computes the value of the stopping criterion*/
  MeshSet *mesh=get_mesh_set();
  bool inner_criterion=inner_criterion_met(*mesh,err);
  
  // does the Feasible Newton iteration until stopping is required.
  while ( nb_iterations<maxIteration && !inner_criterion) {
    
    ++nb_iterations;
  
    // Prints out free vertices coordinates. 
    MSQ_DEBUG_ACTION(3,{
      std::cout << "\n  o Free vertices ("<< num_free_vertices
                <<")original coordinates:\n ";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_CHKERR(err);
      MsqFreeVertexIndexIterator ind1(&pd, err); MSQ_CHKERR(err);
      ind1.reset();
      while (ind1.next()) {
        std::cout << "\t\t\t" << toto1[ind1.value()];
      }
    });
      
    double original_value = 0.0;
    fn_bool=objFunc->evaluate(pd, original_value, err);  MSQ_CHKERR(err);
    if(!fn_bool){
      err.set_msg("Feasible Newton passed invalid patch");
    }
    MSQ_DEBUG_ACTION(3,{std::cout << "  o  original_value: " << original_value
                                  << std::endl;});
    

    // Prints out free vertices coordinates. 
    MSQ_DEBUG_ACTION(3,{
      std::cout << "  o Free vertices new coordinates: \n";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_CHKERR(err);
      MsqFreeVertexIndexIterator ind(&pd, err); MSQ_CHKERR(err);
      ind.reset();
      while (ind.next()) {
        std::cout << "\t\t\t" << toto1[ind.value()];
      }
    });

    bool inner_criterion=inner_criterion_met(*mesh,err);
  }

}


#undef __FUNC__
#define __FUNC__ "FeasibleNewton::terminate_mesh_iteration" 
void FeasibleNewton::terminate_mesh_iteration(PatchData &pd, MsqError &err)
{
  //  std::cout << "- Executing FeasibleNewton::iteration_complete()\n";
}
  
#undef __FUNC__
#define __FUNC__ "FeasibleNewton::cleanup" 
void FeasibleNewton::cleanup()
{
  //  std::cout << "- Executing FeasibleNewton::iteration_end()\n";
}
  

