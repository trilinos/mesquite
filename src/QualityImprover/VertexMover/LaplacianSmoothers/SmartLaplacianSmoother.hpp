/*!
  \file   SmartLaplacianSmoother.hpp
  \brief  

  The SmartLaplacianSmoother Class is the concrete class
  that performs Laplacian Smoothing

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_SmartLaplacianSmoother_hpp 
#define Mesquite_SmartLaplacianSmoother_hpp


#include "Mesquite.hpp"
#include "VertexMover.hpp"

namespace Mesquite
{

  /*   */ 
  class SmartLaplacianSmoother : public VertexMover 
  {
  public:
    SmartLaplacianSmoother(const char* quality_measure);

  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void optimize_nodes_position(PatchData &, 
                                         MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

    char qualityMeasure[80];
  };


}

#endif
