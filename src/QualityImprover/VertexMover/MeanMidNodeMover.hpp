/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#ifndef MESQUITE_MEAN_MID_NODE_MOVER_HPP
#define MESQUITE_MEAN_MID_NODE_MOVER_HPP

/**\file MeanMidNodeMover.hpp
 *\author Jason Kraftcheck
 *\date 2004-12-6
 */

#include "Mesquite.hpp"
#include "Instruction.hpp"

namespace MESQUITE_NS {

class PatchData;

/**\brief Class to adjust positions of higher-order nodes.
 *
 *Move all higher-order nodes to average position of adjacent nodes.
 */
class MeanMidNodeMover : public Instruction
{
public:

  MeanMidNodeMover();
  
  virtual ~MeanMidNodeMover();

  virtual double loop_over_mesh( Mesh* mesh, 
                                 MeshDomain* domain,
                                 MsqError &err);

  virtual std::string get_name() const;

private:

  void fix_mid_nodes( PatchData& pd, MsqError& err );
  
};

} // namespace Mesquite

#endif
