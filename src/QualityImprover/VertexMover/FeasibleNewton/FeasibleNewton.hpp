// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 15-Jan-03 at 08:05:56
//  LAST-MOD:  1-May-03 at 11:14:44 by Thomas Leurent
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
  \author Todd Munson
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

      \brief High Performance implementation of the Feasible Newton algorythm.

      Consider our non-linear objective function
      \f$ f: I\!\!R^{3N} \rightarrow I\!\!R \f$ where \f$ N \f$
      is the number of vertices of the mesh, and \f$ 3N \f$ is therefore the number
      of degrees of freedom of the mesh.
      The Taylor expansion of \f$ f \f$ around the point \f$ x_0 \f$ is 
      \f[ f(x_0+d) = f(x_0) + \nabla f(x_0)d + \frac{1}{2} d^T\nabla^2 f(x_0)d
          + ...  \;\;\; .\f]

      Each iteration of the Newton algorithm tries to find a descent vector that
      minimizes the above quadratic approximation, i.e. it looks for
      \f[ \min_{d} q(d;x_0) = f(x_0) + \nabla f(x_0)d + \frac{1}{2} d^T\nabla^2 f(x_0)d
          \;\; . \f]
      We know that if a quadratic function has a finite minimum, it is reached at the
      point where the function gradient is null and that the function Hessian
      is then positive definite. 
      Therefore we are looking for \f$ d \f$ such that \f$ \nabla q(d;x_0) =0 \f$. We have
      \f[ \nabla q(d;x_0) = \nabla f(x_0) + \nabla^2 f(x_0)d \;\;, \f]
      therefore we must solve for \f$ d \f$ the system
      \f[ \nabla^2 f(x_0)d = -\nabla f(x_0) \;\; . \f]

      We assume that the Hessian is positive definite and we use the conjugate gradient
      algebraic solver to solve the above system. If the conjugate gradient solver finds
      a direction of negative curvature, the Hessian was not positive definite and we take
      a step in that direction of negative curvature, which is a descent direction. 
  */ 
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
