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
    ConjugateGradient(ObjectiveFunction* objective);
      /*!Set the maximum number of search directions to try before moving
        to the next patch.
      */
    void set_maximum_iteration(int iter){
        maxIteration=iter;}

     /*!Set the minimum step length to be allowed before moving to the
       next patch.
      */
    void set_step_size_bound(double step){
      stepBound=step;}

      /*!Set the minimum value of the inifinity norm of the ObjectiveFunction
        gradient to be allowed before moving to the next patch.
      */
    void set_gradient_bound(double grad)
      {
        normGradientBound=grad;
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
    ObjectiveFunction* objFunc;
    int maxIteration;
    double stepBound;
    double normGradientBound;
    Vector3D* fGrad;
    Vector3D* pGrad;
    Vector3D* mCoord;
    Vector3D* fNewGrad;
    int arraySize;
  };

  

}
 
#endif
