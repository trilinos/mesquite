/*!
  \file   ConjugateGradient.hpp
  \brief 

  Conjugate Gradient minimization method ...

  \author Michael Brewer
  \date   2002-06/19
*/

#ifndef Mesquite_ConjugateGradient_hpp 
#define Mesquite_ConjugateGradient_hpp
#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "ObjectiveFunction.hpp"
#include "MsqMeshEntity.hpp"
#include "PatchData.hpp"
#include "QualityImprover.hpp"


namespace Mesquite
{

  /*! \class ConjugateGradient
    \brief Optimizes the objective function using the Polack-Ribiere scheme.
   */ 
  class ConjugateGradient : public VertexMover 
  {
  public:
    ConjugateGradient(ObjectiveFunction* objective, MsqError &err);

    virtual ~ConjugateGradient();
    
      //!Set the patch type
    virtual void set_patch_type(PatchData::PatchType type, MsqError &err, 
                                int patch_param1=0, int patch_param2=0);
      //!Just for debugging purposes or for obtaining more data
      //! during the optimization process.
    void set_debugging_level(int new_lev)
      {
        conjGradDebug=new_lev;
      }
    
  protected:
      
      //!Initialize data for smoothing process
    virtual void initialize(PatchData &pd, MsqError &err);
 
    virtual void optimize_vertex_positions(PatchData &pd, MsqError &err);

    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    
      //!Delete arrays initially created in initialize().
    virtual void cleanup();
   
      //!Returns the step distance to take in the search direction.
    double get_step(PatchData &pd, double f0,int &j, MsqError &err);
      
      //!Culls the vertex list free_vertex_list.     
      //void cull_list(PatchData &pd, double beta, MsqError &err);
    
     
private:
    Vector3D* fGrad;
    Vector3D* pGrad;
    PatchDataVerticesMemento* pMemento;
    Vector3D* fNewGrad;
    int arraySize;
      //just for debugging
    int conjGradDebug;
  };

  

}
 
#endif
