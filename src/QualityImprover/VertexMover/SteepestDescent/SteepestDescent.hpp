/*!
  \file   SteepestDescent.hpp
  \brief  

  The SteepestDescent Class implements the steepest descent algorythm in
  order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QaulityMetric object.

  \author Thomas Leurent
  \date   2002-06-13
*/

#ifndef Mesquite_SteepestDescent_hpp 
#define Mesquite_SteepestDescent_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "ObjectiveFunction.hpp"


namespace Mesquite
{

  /*   */ 
  class SteepestDescent : public VertexMover 
  {
  public:
    SteepestDescent(ObjectiveFunction* of);

    void set_maximum_iteration(int iter){
        maxIteration=iter;}

    void set_lower_gradient_bound(double gradc){
        gradientLessThan=gradc;}
    
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

  private:
    ObjectiveFunction* objFunc;
    double gradientLessThan;
    int maxIteration;
  };
  
}

#endif
