/*!
  \file   SmartLaplacianSmoother.hpp
  \brief  

  The SmartLaplacianSmoother Class implements the SmartLaplacian smoothing
  for a patch with one free vertex. 

  \author Michael Brewer
  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_SmartLaplacianSmoother_hpp 
#define Mesquite_SmartLaplacianSmoother_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "ObjectiveFunction.hpp"
#include <vector>
namespace Mesquite
{

  /*! \class SmartLaplacianSmoother
    Moves free center vertex to the average of the neighboring vertices if
    that move does not decrease the given objective function value.
   */  
  class SmartLaplacianSmoother : public VertexMover 
  {
  public:
    SmartLaplacianSmoother(ObjectiveFunction *obj_func, MsqError &err);
    ~SmartLaplacianSmoother();
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

  private:
    QualityMetric* edgeQM;

    ObjectiveFunction* defaultObjFunc;

  };

  

  
}

#endif
