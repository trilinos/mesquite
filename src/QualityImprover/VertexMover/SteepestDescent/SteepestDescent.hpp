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

  /*! \class SteepestDescent

      This is a very basic implementation of the steepest descent optimization algorythm.
      It works on patches of any size but the step size is hard-wired.
      Obvisouly, this is for testing purposed only. */ 
  class SteepestDescent : public VertexMover 
  {
  public:
    SteepestDescent(ObjectiveFunction* of);


    //! Sets the patch data to global or local.
    virtual void set_patch_type(MeshSet::PatchType type, MsqError &err);

    /*! sets the maximum number of iteration of the steepest descent algorythm,
      i.e. the number of times we compute the gradient and try to move the nodes in the
      opposite direction. This is different from the number of passes over the mesh. */
    void set_maximum_iteration(int iter){
        maxIteration=iter;}

    /*! Sets a minimum value for the gradient. If the gradient is below that value,
      we stop iterating. */  
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
