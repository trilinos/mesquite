/*! \file ShapeImprover.hpp

Header file for the Mesquite::ShapeImprover class

  \author Darryl Melander
  \date   2002-11-14
 */


#ifndef SHAPEIMPROVER_HPP
#define SHAPEIMPROVER_HPP

#include "InstructionWrapper.hpp"

namespace Mesquite
{
  class MeshSet;
  
  /*! \class ShapeImprover
    \brief This class provides a simple interface to enable
           a mesh's shape to be improved.
  */
  class ShapeImprover : public InstructionWrapper
  {
  public:
    ShapeImprover();
    ~ShapeImprover();
    
    void improve_quality(MeshSet &mesh_set, MsqError &err);
    
  private:
    
  };
}

#endif
