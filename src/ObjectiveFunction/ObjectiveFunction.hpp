// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ObjectiveFunction.hpp

Header file for the Mesquite::ObjectiveFunction class

  \author Michael Brewer
  \date   2002-05-23
 */


#ifndef ObjectiveFunction_hpp
#define ObjectiveFunction_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityMetric.hpp"
#include "PatchData.hpp"
#include <list>
#include "MsqMessage.hpp"
namespace Mesquite
{
   class PatchData;
     /*! \class ObjectiveFunction
       \brief Base class for concrete Objective Functions
       ObjectiveFunction contains a pointer to a QualityMetric.  If
       the ObjectiveFunction is associated with more than one
       QualityMetric (i.e., the Objective is a composite, and the
       composed ObjectiveFunctions are associated with different
       QualityMetrics), then the QualityMetric pointer is set
       to NULL..
      */
  class ObjectiveFunction
  {
   public:
      //ObjectiveFunction (QualityMetric *);
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~ObjectiveFunction()
       {};
    
      /*! 
        Evaluate the objective function on a given patch.
      */
    virtual double concrete_evaluate(PatchData &patch, MsqError &err)=0;

      /*! 
        Returns negateFlag*concrete_evaluate(patch, err).
      */
    double evaluate(PatchData &patch, MsqError &err)
       {
         return negateFlag*concrete_evaluate(patch,err);
       }
    
    enum GRADIENT_TYPE{
       NUMERICAL_GRADIENT,
       ANALYTICAL_GRADIENT
    };

      /*!Set gradType to either NUMERICAL_GRADIENT or
        ANALYTICAL_GRADIENT.
      */
    void set_gradient_type(GRADIENT_TYPE grad)
       {
         gradType=grad;
       }

      /*!Calls either compute_numerical_gradient or compute_analytical_gradeint
        depending on the value of gradType.
      */
    void compute_gradient(PatchData &patch, Vector3D *const &grad, MsqError &err, int array_size=-1);          
              
  /*! 
        Return the quality metric associated with this objective function.
        Returns null for composite functions which have multiple associated
        quality metrics.  Use get_quality_metric_list() to retrieve all
        metrics.
      */
    QualityMetric*  get_quality_metric(){
         return qMetric;
       }
    
      /*! 
        returns a list of all associated metrics;
      */
    virtual std::list<QualityMetric*> get_quality_metric_list()
       {
         std::list<QualityMetric*> temp_list;
         temp_list.push_front(qMetric);
         return temp_list;
       }

      //!Set the value of qMetric.
    void set_quality_metric(QualityMetric* qm)
       {
         qMetric=qm;
       }
    
      /*! 
        return the (integer) feasible constraint (1 implies that the
        element is required to stay in the feasible region).
      */
    int get_feasible_constraint(){
         return feasibleFlag;
       }
    
      /*! 
        \brief Set the value of ObjectiveFunction's negateFlag.  Unless
        composite, concrete ObjectiveFunctions should set this flag to
        to the value of the associated QualityMetric's negateFLag.
      */
    void set_negate_flag(int neg)
       {
         negateFlag=neg;
       }

      //!Returns negateFlag
    int get_negate_flag()
       {
         return negateFlag;
       }
    
   protected:

      //!Non-virtual function which numerically computes the gradient
      //!  of the Objective Function.
    void compute_numerical_gradient(PatchData &patch, Vector3D *const &grad,
                                    MsqError &err, int array_size);

     /*! 
        Fills an array of Vector3D, grad, with the gradient of
        the objective function computed using the gradient of the
        quality metric.  If the function has not been over-riden
        in the concrete Objective Function, the base class implementation
        prints a warning and then defaults to numerical gradient.
      */
    virtual void  compute_analytical_gradient(PatchData &patch,
                                              Vector3D *const &grad,
                                              MsqError &err, int array_size){
      PRINT_WARNING("Analytic gradient not implemented for this Objective ",
                    "Function. Defaulting to numerical gradient.\n");
      set_gradient_type(NUMERICAL_GRADIENT);
      compute_numerical_gradient(patch, grad, err, array_size);
    }
    
      //!Returns eps used in the numerical gradient calculation.
    double get_eps(PatchData &pd,int k,MsqVertex* vertex, MsqError &err);
    
      //!Checks to make sure the mesh is in the feasible region.
    int check_feasible(PatchData &pd, MsqError &err);

      //!Sets feasibleFlag.
    void set_feasible(int fflag)
       {
         feasibleFlag=fflag;
       }
    
 private:
    
      
    enum GRADIENT_TYPE gradType;//!Flag for numerical or analytical gradient.
      
    QualityMetric* qMetric;//!Pointer to associated QualityMetric.
      
    int feasibleFlag;//!Set to one if metric has a feasibility constraint.
     
    int negateFlag; /*!Equals one if ObjectiveFunction needs to
        be minimized; equals negative one if ObjectiveFunction needs
        to be maximized.*/
  };

//BEGIN INLINE

#undef __FUNC__
#define __FUNC__ "ObjectiveFunction::compute_gradient"
     /*!  
       Calls either compute_numerical_gradient or compute_analytical_gradient
       depending on the value of gradType.           
      */
   inline void ObjectiveFunction::compute_gradient(PatchData &patch,
                                                   Vector3D *const &grad,
                                                   MsqError &err,
                                                   int array_size)
   {
     switch(gradType){
       case NUMERICAL_GRADIENT:
          compute_numerical_gradient(patch, grad, err, array_size);MSQ_CHKERR(err);
          break;
       case ANALYTICAL_GRADIENT:
          compute_analytical_gradient(patch, grad, err, array_size);MSQ_CHKERR(err);
          break;
     }
   }

} //namespace


#endif // ObjectiveFunction_hpp


