#ifndef MESQUITE_NULL_IMPROVER_HPP
#define MESQUITE_NULL_IMPROVER_HPP
/*!
  \file   NullImprover.hpp
  \brief  The NullImprover Class is a do-nothing VertexMover.  It just
          loops over the mesh without doing any real work.
          It is used to test functions
          found in VertexMover, such as loop_over_mesh().
          
  \author Darryl Melander
  \date   2002-12-10
*/

#include "VertexMover.hpp"

namespace Mesquite
{
  class NullImprover : public VertexMover
  {
  protected:
    virtual void initialize(PatchData &, MsqError &)
      {}
    virtual void cleanup()
      {}
    virtual void optimize_vertex_positions(PatchData &, 
                                           MsqError &)
      {}
    virtual void initialize_mesh_iteration(PatchData &, 
                                           MsqError &)
      {}
    virtual void terminate_mesh_iteration(PatchData &, 
                                          MsqError &)
      {}
  };
}

#endif
