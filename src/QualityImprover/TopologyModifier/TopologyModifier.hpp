/*!
  \file   TopologyModifier.hpp
  \brief  

  The TopologyModifier Class is the base class for swapping algorythms etc ...

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_TopologyModifier_hpp 
#define Mesquite_TopologyModifier_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityImprover.hpp"

namespace Mesquite
{

  /*   */ 
  class TopologyModifier : public QualityImprover 
  {
  public:
    // virtual destructor ensures use of polymorphism during destruction
    virtual ~TopologyModifier()
      {};
    
    virtual void loop_over_mesh(MeshSet &mesh, MsqError &err);

  protected:
    TopologyModifier();

    virtual void iteration_begin() = 0;
    virtual void optimize_connectivity() = 0;
    virtual void iteration_complete() = 0;
    virtual void iteration_end() = 0;

  };
}

#endif
