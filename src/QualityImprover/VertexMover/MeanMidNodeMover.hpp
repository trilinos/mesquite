#ifndef MESQUITE_MEAN_MID_NODE_MOVER_HPP
#define MESQUITE_MEAN_MID_NODE_MOVER_HPP

/**\file MeanMidNodeMover.hpp
 *\author Jason Kraftcheck
 *\date 2004-12-6
 */

#include "Mesquite.hpp"
#include "VertexMover.hpp"

namespace Mesquite {

/**\brief Class to adjust positions of higher-order nodes.
 *
 *Move all higher-order nodes to average position of adjacent nodes.
 */
class MeanMidNodeMover : public VertexMover
{
public:

  MeanMidNodeMover();
  
  virtual ~MeanMidNodeMover();
  
protected:

  virtual void initialize( PatchData& pd, MsqError& err );
  virtual void cleanup();
  virtual void optimize_vertex_positions( PatchData& pd, MsqError& err );
  virtual void initialize_mesh_iteration( PatchData& pd, MsqError& err ); 
  virtual void terminate_mesh_iteration( PatchData& pd, MsqError &err );
  
};

} // namespace Mesquite

#endif
