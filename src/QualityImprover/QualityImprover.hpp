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
namespace Mesquite
{

  class MeshSet;
  class StoppingCriterion;
  
  /*! \class QualityImprover
    Base class for all quality improvers.

  */ 
  class QualityImprover{
  public:
    
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

    /*! \enum culling_method
      Those are the culling method available to the users.
      Developpers: The values used in that enum are used by a bitset,
      so they have to be 2-based (2,4,8,16,32, ...)
      */
    enum culling_method {
      NO_BOUNDARY_VTX = 1, /*!< removes vertices on the boundary. (i.e. with a TSTT tag "boundary"). */
      CULL_METHOD_2 = 2,   /*!< no other culling method yet. */
      CULL_METHOD_3 = 4,
      CULL_METHOD_4 = 8
    };
    //! Sets on a culling criterion.
    void add_culling_method(enum culling_method cm);
    //! No culling performed (sets off all culling criteria).
    void no_culling_method();
    //! Sets off a certain culling criteria. 
    void remove_culling_method(enum culling_method cm);

    void set_stopping_criterion(StoppingCriterion* crit)
      {
        stoppingCriterion=crit;
      }

    void set_patch_depth(int depth) { patchDepth = depth; }
    int get_patch_depth() { return patchDepth; }
    
  protected:

    QualityImprover() : mMeshSet(0), qualityImproverName("noname"),
                        patchDepth(1), cullingMethodBits(0),
                        stoppingCriterion(0) {}
    
    friend class MeshSet;
      //friend double QualityMetric::evaluate_element(MsqMeshEntity* element, MsqError &err);
      //friend double QualityMetric::evaluate_node(MsqNode* node, MsqError &err);
    long unsigned int get_culling_method_bits() { return cullingMethodBits; }
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
    long unsigned int cullingMethodBits;
    StoppingCriterion* stoppingCriterion;
  };

#undef __FUNC__
#define __FUNC__ "QualityImprover::inner_criterion_met"
  /*! \fn QualityImprover::inner_criterion_met(MeshSet &ms, MsqError &err)
    */
  inline bool QualityImprover::inner_criterion_met(MeshSet &ms, MsqError &err)
  {
    
    bool inner_criterion_met;
    int depth=get_patch_depth();
    StoppingCriterion* crit = get_stopping_criterion();
    if(crit==0){
      err.set_msg("Stopping Criterion pointer is Null");
      return true;
    }
    if(depth<1)
      inner_criterion_met=crit->stop(ms,err);
    else
      inner_criterion_met=0;
    return inner_criterion_met;
  }


#undef __FUNC__
#define __FUNC__ "QualityImprover::add_culling_method"
  /*! \fn QualityImprover::add_culling_method(enum culling_method cm)
    */
  inline void QualityImprover::add_culling_method(enum culling_method cm)
  {
    cullingMethodBits |= cm;
  }
  
#undef __FUNC__
#define __FUNC__ "QualityImprover::no_culling_method"
  /*! \fn QualityImprover::no_culling_method()
   */
  inline void QualityImprover::no_culling_method()
  {
    cullingMethodBits = 0;
  }
  
#undef __FUNC__
#define __FUNC__ "QualityImprover::remove_culling_method"
  /*! \fn QualityImprover::remove_culling_method(enum culling_method cm)
    */
  inline void QualityImprover::remove_culling_method(enum culling_method cm)
  {
    cullingMethodBits &= ~cm;
  }

}

#endif
