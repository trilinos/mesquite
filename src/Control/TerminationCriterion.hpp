// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file TerminationCriterion.hpp

Header file for the TerminationCriterion classes.

  \author Michael Brewer
  \date   Feb. 14, 2003
 */


#ifndef TerminationCriterion_hpp
#define TerminationCriterion_hpp

#include <string>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "PatchData.hpp"
#include "ObjectiveFunction.hpp"
#include "MsqTimer.hpp"

namespace Mesquite {
   class MeshSet;
  /*! \class TerminationCriterion

      \brief The TerminationCriterion class contains
  */
  class TerminationCriterion
  {
  public:
    
     /*! */
    enum TCType {
       NONE    = 0,
       CRIT_1  = 1,
       CRIT_2  = 2, 
       QUALITY_IMPROVEMENT = 4,
       ITERATION_BOUND  = 8,
       CPU_TIME  = 16,
       VERTEX_MOVEMENT  = 32,
       SUCCESSIVE_IMPROVEMENTS  = 64,
       CRIT_8  = 128
    };

      //!Constructor which does not take any arguements
    TerminationCriterion();
    
      //!Destructor
    ~TerminationCriterion(){};

      //Functions with which the user can specify the criteria to be used
      //!Sets the criterion by specifing the TCType and the eps value
    void add_criterion_type_with_double(TCType tc_type, double eps,
                                        MsqError &err);
      //!Sets the criterion by specifing the TCType and the integer value
    void add_criterion_type_with_int(TCType tc_type, int bound,
                                     MsqError &err);
      //!Sets the type of criterion that the user would like to
      //! use for culling purposes (along with the associated tolerance.
    void set_culling_type(TCType tc_type, double eps, MsqError &err);
    
      //Functions usually called from vertex mover (either concrete or base)
      //!Does preliminary set up for the TerminationCriterion
    void initialize(MeshSet &ms, PatchData &pd, MsqError &err);
      //!Does some setup and initialization so the criterion can be used
      //! (not necessarily for the first time).
    void reset(PatchData &pd, ObjectiveFunction* obj_ptr, MsqError &err);
      //!reset called using a MeshSet object instead of PatchData.
    void reset(MeshSet &ms, ObjectiveFunction* obj_ptr, MsqError &err);
      
     //!Returns true if termination criterion is met (for the inner loop). 
    bool terminate(PatchData &pd, ObjectiveFunction* obj_ptr, MsqError &err);
      //!Returns true if termination criterion is met (for the outer loop).
    bool terminate(MeshSet &ms, ObjectiveFunction* obj_ptr, MsqError &err);
      //!Function which determines whether this patch should be 'culled'
    bool cull_vertices(PatchData &pd, ObjectiveFunction* obj_ptr,
                       MsqError &err);
      //!Cleans up after the TerminationCriterion is finished.
    void cleanup(MeshSet &ms, MsqError &err);
      
    
 protected:
    
 private:
    long unsigned int terminationCriterionFlag;//!<Bit flag of termination crit
    long unsigned int cullingMethodFlag;/*!<Bit flag of criterion for culling*/
    long unsigned int totalFlag;/*!<Bit flag for both culling and terminating.*/
      //epsiloon used in culling methods.
    double cullingEps;

      //Data not specific to a single criterion
    double initialOFValue;
    double previousOFValue;
    double lowerOFBound;
      //if we need to create a global patch from a meshset.
    PatchDataParameters globalPatchParams;
      //Data specific to termination criterion 1 (gradient bounds)
    Vector3D* mGrad;
    double initialGradNorm;
    double crit1Eps;
      //Data specific to termination criterion 2 (KKT)
      //???????????????????????????????????????????
      //Data specific to termination criterion 3 (Quality Improvement)
    double qualityImprovementEps;
      //Data specific to termination criterion 4 (inner iterations)
    int iterationBound;
    int iterationCounter;
      //Data specific to termination criterion 5 (cpu time)
    Timer mTimer;
    double timeBound;
      //Data specific to termination criterion 6 (vertex movement)
    PatchDataVerticesMemento* initialVerticesMemento;
      //PatchDataVerticesMemento* previousVerticesMemento;//if we want relative
    double vertexMovementEps;
      //Data specific to termination criterion 7 (successive improvement to F)
    double successiveImprovementsEps;
      //crit 8
    double crit8Bound;
  };

   
} //namespace


#endif // TerminationCriterion_hpp
