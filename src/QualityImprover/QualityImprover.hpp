/*!
  \file   QualityImprover.hpp
  \brief  

  The Quality Improver Class is the base class for all the algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_QualityImprover_hpp 
#define Mesquite_QualityImprover_hpp

#include <string>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "StoppingCriterion.hpp"
#include "PatchData.hpp"
#include "MeshSet.hpp"

namespace Mesquite
{

//   class MeshSet;
//   class StoppingCriterion;
  
  /*! \class QualityImprover
    Base class for all quality improvers.

  */ 
  class QualityImprover : public PatchDataUser {
  public:

    // Constructor is protected ... see below.
    
     // virtual destructor ensures use of polymorphism during destruction
    virtual ~QualityImprover() { };
    
    virtual void loop_over_mesh(MeshSet &ms, MsqError &err) = 0;

    //! provides a name to the QualityImprover (use it in constructor).
    void set_name(std::string name)
      {
        qualityImproverName = name;
      };
    
    //! retrieves the QualityImprover name. A default name should be set in the constructor.
    std::string get_name() { return qualityImproverName; }

    void set_stopping_criterion(StoppingCriterion* crit)
      {
        stoppingCriterion=crit;
      }

  protected:

    /*! The default constructor initialises a few member variables
        to default values.
        This can be reused by concrete class constructor. */
    QualityImprover() : mMeshSet(0), qualityImproverName("noname"),
                        stoppingCriterion(0) {}
    
    friend class MeshSet;
      //friend double QualityMetric::evaluate_element(MsqMeshEntity* element, MsqError &err);
      //friend double QualityMetric::evaluate_node(MsqNode* node, MsqError &err);
    StoppingCriterion* get_stopping_criterion() { return stoppingCriterion; }
    bool inner_criterion_met(MeshSet &ms, MsqError &err);
    
    const MeshSet* get_mesh_set() const
      { return mMeshSet; }
    MeshSet* get_mesh_set()
      { return mMeshSet; }
    
    
    void set_mesh_set(MeshSet *ms)
      {
        mMeshSet=ms;
      }
    
  private:
    MeshSet* mMeshSet;
    std::string qualityImproverName;
    int patchDepth;
    StoppingCriterion* stoppingCriterion;
  };

#undef __FUNC__
#define __FUNC__ "QualityImprover::inner_criterion_met"
  /*! \fn QualityImprover::inner_criterion_met(MeshSet &ms, MsqError &err)
    */
  inline bool QualityImprover::inner_criterion_met(MeshSet &ms, MsqError &err)
  {
    
    bool inner_criterion_met;
    StoppingCriterion* crit = get_stopping_criterion();
    if(crit==0){
      err.set_msg("Stopping Criterion pointer is Null");
      return true;
    }
    if(this->get_patch_type()==PatchData::GLOBAL_PATCH)
      inner_criterion_met=crit->stop(ms,err);
    else
      inner_criterion_met=0;
    return inner_criterion_met;
  }

}

#endif
