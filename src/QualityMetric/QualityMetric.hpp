// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file QualityMetric.hpp
    \brief
Header file for the Mesquite::QualityMetric class

  \author Thomas Leurent
  \date   2002-05-01
 */


#ifndef QualityMetric_hpp
#define QualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include <math.h>
#include <string.h>

namespace Mesquite
{
   
     /*! \class QualityMetric
       \brief Base class for concrete quality metrics.
     */
   class MsqVertex;
   class MsqMeshEntity;
   class PatchData;
   
   class QualityMetric
   {
   public:
       // This is defined in each concrete class.  It isn't virtual, so
       // it doesn't exist in the base class.
       //   static void create_new() = 0;
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~QualityMetric()
        {};
     
       /*!EvaluationMode allows you to indicate whether we are evaluating
         the metric based on vertices, element vertices, or element gauss
         points.
       */
     enum EvaluationMode
     {
       VERTEX,
       ELEMENT_VERTICES,
       LINEAR_GAUSS_POINTS,
       QUADRATIC_GAUSS_POINTS,
       CUBIC_GAUSS_POINTS
     };
       //!Set the evuation mode for the metric
     void set_mode(EvaluationMode mode, MsqError &err)
        {
          switch(mode)
          {
            case(VERTEX):
            case(ELEMENT_VERTICES):
            case(LINEAR_GAUSS_POINTS):
            case(QUADRATIC_GAUSS_POINTS):
            case(CUBIC_GAUSS_POINTS):
               evalMode=mode;
               break;
            default:
               err.set_msg("Requested Mode (sample point locations)  Not Implemented");
          }
          return;
        }
     
       //!Returns the evaluation mode for the metric
     inline EvaluationMode  get_evaluation_mode()
        { return evalMode; }
     
       /*!AveragingMethod allows you to set how the quality metric values
         attained at each sample point will be averaged together to produce
         a single metric value for an element.
       */
     enum AveragingMethod
     {
        NONE,
        LINEAR,
        RMS,
        MINIMUM,
        MAXIMUM,
        HARMONIC,
        GEOMETRIC,
        SUM
     };
     
       /*!Set the averaging method for the quality metric. Current
         options are
         NONE: the values are not averaged,
         GEOMETRIC:  the geometric average,
         HARMONIC:  the harmonic average,
         LINEAR:  the linear average,
         MAXIMUM:  the maximum value,
         MINIMUM:  the minimum value,
         RMS:  the root-mean-squared average,
         SUM:  the sum of the values
       */
     inline void set_averaging_method(AveragingMethod method, MsqError &err)
        {
          switch(method)
          {
            case(NONE):
            case(GEOMETRIC):
            case(HARMONIC):
            case(LINEAR):
            case(MAXIMUM):
            case(MINIMUM):
            case(RMS):
            case(SUM):
               avgMethod=method;
               break;
            default:
               err.set_msg("Requested Averaging Method Not Implemented");
          };
          return;
        }
     
       /*! Set feasible flag (i.e., does this metric have a feasible region
         that the mesh must maintain.)
       */
     inline void set_feasible_constraint(int alpha)
        { feasible=alpha; }
     
       //!Returns the feasible flag for this metric
     inline int get_feasible_constraint()
        { return feasible; }
     
       //!Sets the name of this metric
     inline void set_name(std::string st)
        { metricName=st; }
     
       //!Returns the name of this metric (as a string).
     inline std::string get_name()
        { return metricName; }
     
       //!Evaluate the metric for a vertex
     virtual double evaluate_vertex(PatchData &pd, MsqVertex* vertex,
                                    MsqError &err)
        {
          err.set_msg("No implementation for a "
                      "vertex-version of this metric.");
          return 0.0;
        }
     
       //!Evaluate the metric for an element
     virtual double evaluate_element(PatchData &pd,
                                     MsqMeshEntity* element,
                                     MsqError &err)
        {
          err.set_msg("No implementation for a element-version of this metric.");
          return 0.0;
        }
     
       /*!\enum GRADIENT_TYPE Sets to either NUMERICAL_GRADIENT or
         ANALYTICAL_GRADINET*/
     enum GRADIENT_TYPE
     {
        NUMERICAL_GRADIENT,
        ANALYTICAL_GRADIENT
     };
     
       //!Sets gradType for this metric.
     void set_gradient_type(GRADIENT_TYPE grad)
        { gradType=grad; }
     
       /*!Calls either compute_numerical_gradient or
         compute_analytical_gradient for gradType equal
         NUMERICAL_GRADIENT or ANALYTICAL_GRADIENT, respectively.
       */
     void compute_gradient(PatchData &pd, MsqMeshEntity* element,
                           MsqVertex &vertex, Vector3D &grad_vec,
                           MsqError &err);
     
       /*! Set the value of QualityMetric's negateFlag.  Concrete
         QualityMetrics should set this flag to -1 if the QualityMetric
         needs to be maximized.
       */
     void set_negate_flag(int neg)
        { negateFlag=neg; }
     
       //!Returns negateFlag.
     int get_negate_flag()
        { return negateFlag; }
     
  protected:
     
       /*!Constructor defaults concrete QualityMetric's to
         gradType=NUMERCIAL_GRADIENT and negateFlag=1.
         Concrete QualityMetric constructors over-write these defaults
         when appropriate.
       */
     QualityMetric() :
         gradType(NUMERICAL_GRADIENT),
         negateFlag(1)
        {}
     
       /*! average_metrics takes an array of length num_values and averages the
         contents using averaging method 'avgMethod'.
       */
     double average_metrics(double* metric_values, int num_values,
                            MsqError &err);
     
       /*!Non-virtual function which numerically computes the gradient
         of a QualityMetric of a given element for a given free vertex.
       */
     void compute_numerical_gradient(PatchData &pd, MsqMeshEntity* element,
                                     MsqVertex &vertex, Vector3D &grad_vec,
                                     MsqError &err);
     
       /*Virtual function that computes the gradient of the QualityMetric
         analytically.  The base class implementation of this function
         simply prints a warning and calls compute_numerical_gradient
         to calculate the gradient. 
       */
     virtual void compute_analytical_gradient(PatchData &pd,
                                              MsqMeshEntity* element,
                                              MsqVertex &vertex,
                                              Vector3D &grad_vec,
                                              MsqError &err);
     
     friend class MsqMeshEntity;

     // TODO : pass this private and write protected access fucntions.
     AveragingMethod avgMethod;
     EvaluationMode evalMode;
     int feasible;
     std::string metricName;
  private:
     GRADIENT_TYPE gradType;
     int negateFlag;
   };

#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_gradient"
/*! 
  \brief Calls compute_numerical_gradient if gradType equals
  NUMERCIAL_GRADIENT.  Calls compute_analytical_gradient if 
  gradType equals ANALYTICAL_GRADIENT;
*/
   inline void QualityMetric::compute_gradient(PatchData &pd,
                                               MsqMeshEntity* element,
                                               MsqVertex &vertex,
                                               Vector3D &grad_vec,
                                               MsqError &err)
   {
     switch(gradType)
     {
       case NUMERICAL_GRADIENT:
          compute_numerical_gradient(pd, element, vertex, grad_vec, err);
          MSQ_CHKERR(err);
          break;
       case ANALYTICAL_GRADIENT:
          compute_analytical_gradient(pd, element, vertex, grad_vec, err);
          MSQ_CHKERR(err);
          break;
     }
   }
   
   
#undef __FUNC__
#define __FUNC__ "QualityMetric::average_metrics"
     /*! 
       average_metrics takes an array of length num_value and averages the
       contents using averaging method 'method'.
     */
   inline double QualityMetric::average_metrics(double* metric_values,
                                                int num_values, MsqError &err)
   {
     double MSQ_MAX=1e10;
     double total_value=0.0;
     int i=0;
       //if no values, return zero
     if (num_values==0){
       return 0.0;
     }
     
     switch(avgMethod){
       case GEOMETRIC:
          total_value=1.0;
          for (i=0;i<num_values;++i){
            total_value*=(metric_values[i]);
          }
          total_value=pow(total_value, (1/((double) num_values)));
          break;
          
       case HARMONIC:
            //ensure no divide by zero, return zero
          for (i=0;i<num_values;++i){
            if(metric_values[i]<MSQ_MIN){
              if(metric_values[i]>MSQ_MIN){
                return 0.0;
              }
            }
            total_value+=(1/metric_values[i]);
          }
            //ensure no divide by zero, return MSQ_MAX_CAP
          if(total_value<MSQ_MIN){
            if(total_value>MSQ_MIN){
              return MSQ_MAX_CAP;;
            }
          }
          total_value=num_values/total_value;
          break;
          
       case LINEAR:
          for (i=0;i<num_values;++i){
            total_value+=metric_values[i];
          }
          total_value/=num_values;
          break;
          
       case MAXIMUM:
          for (i=0;i<num_values;++i){
            if(metric_values[i]>total_value){
              total_value=metric_values[i];
            }
          }
          break;
          
       case MINIMUM:
          total_value=MSQ_MAX;
          for (i=0;i<num_values;++i){
            if(metric_values[i]<total_value){
              total_value=metric_values[i];
            }
          }
          break;
          
       case NONE:
          err.set_msg("Averaging method set to NONE");
          break;
          
       case RMS:
          for (i=0;i<num_values;++i){
            total_value+=(metric_values[i]*metric_values[i]);
          }
          total_value/=num_values;
          total_value=sqrt(total_value);
          break;
          
       case SUM:
          for (i=0;i<num_values;++i){
            total_value+=metric_values[i];
          }
          break;
          
       default:
            //Return error saying Averaging Method mode not implemented
          err.set_msg("Requested Averaging Method Not Implemented");
          return 0;
     }
     return total_value;
   }
   
   
} //namespace


#endif // QualityMetric_hpp
