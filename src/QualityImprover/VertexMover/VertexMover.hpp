/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   VertexMover.hpp
  \brief  

  The VertexMover Class is the base class for all the smoothing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_VertexMover_hpp 
#define Mesquite_VertexMover_hpp


#include "Mesquite.hpp"
#include "QualityImprover.hpp"
#include "OFEvaluator.hpp"

namespace MESQUITE_NS
{
  class ObjectiveFunction;
  class PatchData;
  class ParallelMesh;

  /*! \class VertexMover
    Base class for all Vertex Movers.
   */  
  class VertexMover : public QualityImprover 
  {
  protected:
    VertexMover( ObjectiveFunction* OF = NULL, bool Nash = true );
    
  public:
    // virtual destructor ensures use of polymorphism during destruction
    virtual ~VertexMover() { };
    
    virtual double loop_over_mesh( Mesh* mesh, 
                                   MeshDomain* domain,
                                   const Settings* settings,
                                   MsqError &err);

    virtual double loop_over_mesh( ParallelMesh* mesh, 
                                   MeshDomain* domain,
                                   const Settings* settings,
                                   MsqError &err);

  protected:

    virtual void initialize(PatchData &pd, MsqError &err) = 0;
    virtual void cleanup() = 0;
    virtual void optimize_vertex_positions(PatchData &pd, 
                                           MsqError &err) = 0; // modifies the PatchData object

    virtual void initialize_mesh_iteration(PatchData &pd, 
                                         MsqError &err) = 0;
    virtual void terminate_mesh_iteration(PatchData &, 
                                         MsqError &err) = 0;

  
    OFEvaluator& get_objective_function_evaluator()
      { return objFuncEval; }
  
  private:
    OFEvaluator objFuncEval;
  };

  
} // namespace
#endif // Mesquite_VertexMover_hpp
