// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 15-Jan-03 at 08:05:56
//  LAST-MOD: 20-Feb-03 at 13:18:31 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*!
  \file   FeasibleNewton.hpp
  \brief  

  The FeasibleNewton Class implements the newton non-linear programming algorythm
  in order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QualityMetric object.

  \author Thomas Leurent
  \date   2003-01-15
*/
// DESCRIP-END.
//

#ifndef MSQ_FeasibleNewton_hpp 
#define MSQ_FeasibleNewton_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "ObjectiveFunction.hpp"


namespace Mesquite
{

  /*! \class FeasibleNewton

      High Performance algorythm of the Feasible Newton algorythm (tmunson@mcs.anl.gov). */ 
  class FeasibleNewton : public VertexMover 
  {
  public:
    FeasibleNewton(ObjectiveFunction* of);


    /*! sets the maximum number of iteration of the Feasible Newton algorythm. */
    void set_maximum_iteration(int iter){
      maxIteration=iter;}

    /*! Sets a minimum value for the gradient. If the gradient is below that value,
      we stop iterating. */  
    void set_lower_gradient_bound(double gradc){
        convTol=gradc;}
    
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                           MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

  private:
    double convTol;
    int maxIteration;
    MsqHessian mHessian;
    PatchDataVerticesMemento* coordsMem;
  };
  
}

#endif // MSQ_FeasibleNewton_hpp 
