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
       GRADIENT_NORM_ABSOLUTE = 1,
       GRADIENT_NORM_RELATIVE = 2,
       KKT  = 4, 
       QUALITY_IMPROVEMENT_ABSOLUTE = 8,
       QUALITY_IMPROVEMENT_RELATIVE = 16,
       NUMBER_OF_ITERATES = 32,
       CPU_TIME  = 64,
       VERTEX_MOVEMENT_ABSOLUTE  = 128,
       VERTEX_MOVEMENT_RELATIVE  = 256,
       SUCCESSIVE_IMPROVEMENTS_ABSOLUTE = 512,
       SUCCESSIVE_IMPROVEMENTS_RELATIVE = 1024,
       BOUNDED_VERTEX_MOVEMENT = 2048
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
      //Private member funcitons
      //!Currently computes the L_inf norm of an array of Vector3D of
      //!length len.
    double compute_gradient_norm(Vector3D * const vec, int len, MsqError &err);

      //PRIVATE DATA MEMBERS
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
    int gradSize;
    double initialGradNorm;
    double gradNormAbsoluteEps;
    double gradNormRelativeEps;
      //Data specific to termination criterion 2 (KKT)
      //???????????????????????????????????????????
      //Data specific to termination criterion 3 (Quality Improvement)
    double qualityImprovementAbsoluteEps;
    double qualityImprovementRelativeEps;
      //Data specific to termination criterion 4 (inner iterations)
    int iterationBound;
    int iterationCounter;
      //Data specific to termination criterion 5 (cpu time)
    Timer mTimer;
    double timeBound;
      //Data specific to termination criterion 6 (vertex movement)
    PatchDataVerticesMemento* initialVerticesMemento;
    PatchDataVerticesMemento* previousVerticesMemento;//if we want relative
    double vertexMovementAbsoluteEps;
    double vertexMovementRelativeEps;
    
      //Data specific to termination criterion 7 (successive improvement to F)
    double successiveImprovementsAbsoluteEps;
    double successiveImprovementsRelativeEps;
      //crit 8
    double boundedVertexMovementEps;
  };

   
} //namespace


#endif // TerminationCriterion_hpp