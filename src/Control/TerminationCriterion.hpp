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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file TerminationCriterion.hpp

Header file for the TerminationCriterion classes.

  \author Michael Brewer
  \author Thomas Leurent
  \date   Feb. 14, 2003
 */


#ifndef TerminationCriterion_hpp
#define TerminationCriterion_hpp

#include "Mesquite.hpp"
#include "PatchDataUser.hpp"
#include "MsqTimer.hpp"

#include <string>

namespace Mesquite
{
   class MeshSet;
   class MsqError;
   class ObjectiveFunction;
   
  /*! \class TerminationCriterion

      \brief The TerminationCriterion class contains functionality to
      terminate the VertexMover's optimization.

      The TerminationCriterion class has three roles.  It
      is used to terminate the optimization on a single patch; it
      is used to terminate the iterations over all patches in the
      mesh; and it is used to cull vertices frm the optimization
      processes.  Thus, for each optimzation, two TerminationCriterion
      objects are used.  The class contains five important member
      functions used in the VertexMover:  initialize(), reset(),
      terminate(), cull_vertices(), and cleanup().  These functions
      are each explained in detail below.  In general, the only one
      of these functions called directly from a concrete VertexMover
      is terminate() which allows the concrete VertexMover to determine
      when to stop producing new iterates on a given patch.  All other
      functionality is handled from the base VertexMover base class.

      There are several different types of termination criteria
      available.  These types are listed in teh enumberation
      TCType.  Multiple criteria types can be set on a given
      TermiantionCriterion object, and when this occurs, the
      optimization process will terminate whenever any of the
      criteria have been satisfied.
      
  */
  class TerminationCriterion
  {
  public:
    
     /*! \enum TCType  defines the termination criterion */
    enum TCType {
       NONE    = 0,
       //! checks the gradient \f$\nabla f \f$ of objective function 
       //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
       //! and stops when \f$\sqrt{\sum_{i=1}^{3N}\nabla f_i^2}<d\f$  
       GRADIENT_L2_NORM_ABSOLUTE = 1,  
       //! checks the gradient \f$\nabla f \f$ of objective function 
       //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
       //! and stops when \f$ \max_{i=1}^{3N} \nabla f_i < d \f$  
       GRADIENT_INF_NORM_ABSOLUTE = 2,
         //!terminates on the j_th iteration when
         //! \f$\sqrt{\sum_{i=1}^{3N}\nabla f_{i,j}^2}<d\sqrt{\sum_{i=1}^{3N}\nabla f_{i,0}^2}\f$
         //!  That is, terminates when the norm of the gradient is small
         //! than some scaling factor times the norm of the original gradient. 
       GRADIENT_L2_NORM_RELATIVE = 4,
       //!terminates on the j_th iteration when
         //! \f$\max_{i=1 \cdots 3N}\nabla f_{i,j}<d \max_{i=1 \cdots 3N}\nabla f_{i,0}\f$
         //!  That is, terminates when the norm of the gradient is small
         //! than some scaling factor times the norm of the original gradient.
         //! (Using the infinity norm.)
       GRADIENT_INF_NORM_RELATIVE = 8,
         //! Not yet implemented.
       KKT  = 16,
         //!Terminates when the objective function value is smaller than
         //! the given scalar value.
       QUALITY_IMPROVEMENT_ABSOLUTE = 32,
         //!Terminates when the objective function value is smaller than
         //! the given scalar value times the original objective function
         //! value.
       QUALITY_IMPROVEMENT_RELATIVE = 64,
         //!Terminates when the number of iterations exceeds a given integer.
       NUMBER_OF_ITERATES = 128,
         //!Terminates when the algorithm exceeds an allotted time limit
         //! (given in seconds).
       CPU_TIME  = 256,
         //!Terminates when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value.
       VERTEX_MOVEMENT_ABSOLUTE  = 512,
         //!Terminates when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value
         //! times the maximum distance moved by any vertex over the
         //! entire course of the optimization.
       VERTEX_MOVEMENT_RELATIVE  = 1024,
         //!Terminates when the decrease in the objective function value since
         //! the previous iteration is below the given value.
       SUCCESSIVE_IMPROVEMENTS_ABSOLUTE = 2048,
         //!Terminates when the decrease in the objective function value since
         //! the previous iteration is below the given value times the
         //! decrease in the objective function value since the beginning
         //! of this optimization process.
       SUCCESSIVE_IMPROVEMENTS_RELATIVE = 4096,
         //!Terminates when any vertex leaves the bounding box, defined
         //! by the given value, d.  That is, when the absolute value of
         //! a single coordinate of vertex's position exceeds d.
       BOUNDED_VERTEX_MOVEMENT = 8192
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
      //!Removes the criterion by specifing just the TCType.
    void remove_criterion_type(TCType tc_type, MsqError &err);
    
      //!Sets the type of criterion that the user would like to
      //! use for culling purposes (along with the associated tolerance.
    void set_culling_type(TCType tc_type, double eps, MsqError &err);
      //!Removes any previously set culling types (sets the culling
      //! type to be NONE).
    void remove_culling(MsqError &err);
    
      //Functions usually called from vertex mover (either concrete or base)
      //!Does preliminary set up for the TerminationCriterion
    void initialize(MeshSet &ms, PatchData &pd, MsqError &err);
      //!Does some setup and initialization so the criterion can be used
      //! (not necessarily for the first time).  When it makes sense,
      //! the termination criteria are checked.  The return value
      //! is similar to that of terminate().
    bool reset(PatchData &pd, ObjectiveFunction* obj_ptr, MsqError &err);
      //!reset called using a MeshSet object instead of PatchData.
      //! When it makes sense,
      //! the termination criteria are checked. 
    bool reset(MeshSet &ms, ObjectiveFunction* obj_ptr, MsqError &err);
      
     //!Returns true if termination criterion is met (for the inner loop). 
    bool terminate(PatchData &pd, ObjectiveFunction* obj_ptr, MsqError &err);
      //!Returns true if termination criterion is met (for the outer loop).
    bool terminate(MeshSet &ms, ObjectiveFunction* obj_ptr, MsqError &err);
      //!Returns true if termination criterion is met (for the inner loop).
      //!  Also supplies the function and gradient values for effeciency.
    bool terminate_with_function_and_gradient(PatchData &pd,
                                              ObjectiveFunction* obj_ptr,
                                              double func_val,
                                              Vector3D* sup_grad,
                                              MsqError &err);
    
      //!Function which determines whether this patch should be 'culled'
    bool cull_vertices(PatchData &pd, ObjectiveFunction* obj_ptr,
                       MsqError &err);
      //!Cleans up after the TerminationCriterion is finished.
    void cleanup(MeshSet &ms, MsqError &err);

      //!This function returns the current function value.
      /*! \todo Michael:  this function is not reliable.  It
        needs to be more robust.  How do we know whether
        currentOFValue got updated or not?  We may want to
        make sure that all the criteria get checked.*/
    double get_current_function_value()
       {return currentOFValue;}
    
 protected:
    
 private:
    //PRIVATE DATA MEMBERS
    long unsigned int terminationCriterionFlag;//!<Bit flag of termination crit
    long unsigned int cullingMethodFlag;/*!<Bit flag of criterion for culling*/
    long unsigned int totalFlag;/*!<Bit flag for both culling and terminating.*/
      //epsiloon used in culling methods.
    double cullingEps;

      //Data not specific to a single criterion
    double initialOFValue;
    double previousOFValue;
    double currentOFValue;
    double lowerOFBound;
      //if we need to create a global patch from a meshset.
    PatchDataParameters globalPatchParams;
      //Data specific to termination criterion 1 (gradient bounds)
    Vector3D* mGrad;
    int gradSize;
    double initialGradL2Norm;
    double gradL2NormAbsoluteEps;
    double gradL2NormRelativeEps;
    double initialGradInfNorm;
    double gradInfNormAbsoluteEps;
    double gradInfNormRelativeEps;
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

      //Variables for usr supplied data
      //true if function was supplied
    bool functionSupplied;
      //true if function was supplied
    bool gradientSupplied;
      //place holder for supplied function value is stored in currentOFValue
 
      //place holder for supplied Gradient
    Vector3D* suppliedGradientArray;
    
  };

} //namespace


#endif // TerminationCriterion_hpp
