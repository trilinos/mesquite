// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ConditionNumberQualityMetric.hpp

Header file for the Mesquite::ConditionNumberQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef ConditionNumberQualityMetric_hpp
#define ConditionNumberQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ShapeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
//Michael delete
#include "MsqMessage.hpp"

namespace Mesquite
{
     /*! \class ConditionNumberQualityMetric
       \brief Computes the condition number of given element.
       The``condition number" is scaled between one and infinity,
       with an ideal element having condition number one.
     */
   class ConditionNumberQualityMetric : public ShapeQualityMetric
   {
  public:
 
       /*!Returns a poitner to a ShapeQualityMetric.  Defaults the metric
         to the LINEAR averaging method and ELEMENT_VERTICES sample points.
         The default metric name is "Condition Number".  It does require a
         feasible region, and the metric needs to be maximized.
       */
     static ShapeQualityMetric* create_new(){
       
       ShapeQualityMetric* m = new ConditionNumberQualityMetric();
       return m;
     }
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~ConditionNumberQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     double evaluate_element(PatchData &pd, MsqMeshEntity *element,
                             MsqError &err);
     
       // evaluate using patch data raw arrays 
       //double evaluate_element(PatchData *pd, int element_index,
       //                    MsqError &err);

     //! Evaluate the "condition number" for a vertex
     double evaluate_vertex(PatchData &pd, MsqVertex *vertex, MsqError &err);
     
  protected:

     double compute_condition_number(Vector3D* jacobian_vectors,
                                     int num_jacobian_vectors,
                                     MsqError &err);
     
     
  private:
     
     ConditionNumberQualityMetric();
    
  };
//BEGIN INLINE FUNCTIONS

   inline double ConditionNumberQualityMetric::compute_condition_number(
      Vector3D* jacobian_vectors, int num_jacobian_vectors, MsqError &err)
   {
       //PRINT_INFO("INSIDE CONDITION NUMBER COMPUTE_CON\n");
     double temp_var=0;
     double return_val=0;
     if(num_jacobian_vectors==2){
       temp_var=fabs((jacobian_vectors[0]*jacobian_vectors[1]).length());
       return_val=jacobian_vectors[0].length_squared();
       return_val+=jacobian_vectors[1].length_squared();
       if(temp_var>=MSQ_MIN){ //if not degenerate
         return_val/=2*temp_var;
       }
       else{
         return_val=MSQ_MAX_CAP;
       }
       return return_val;
     }
     
       //if three jacobian vectors (3D elem)
     else if(num_jacobian_vectors==3){
         //norm squared of J
       double term1=jacobian_vectors[0]%jacobian_vectors[0]+
          jacobian_vectors[1]%jacobian_vectors[1]+
          jacobian_vectors[2]%jacobian_vectors[2];
         //norm squared of adjoint of J
       double term2=(jacobian_vectors[0]*jacobian_vectors[1])%
          (jacobian_vectors[0]*jacobian_vectors[1])+
          (jacobian_vectors[1]*jacobian_vectors[2])%
          (jacobian_vectors[1]*jacobian_vectors[2])+
          (jacobian_vectors[2]*jacobian_vectors[0])%
          (jacobian_vectors[2]*jacobian_vectors[0]);
         //det of J
       temp_var=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
       return_val=sqrt(term1*term2);
       if(return_val>MSQ_MIN){
           //if not degenerate or inverted???
         return_val/=(3*temp_var);
       }
       else{
         return_val=MSQ_MAX_CAP;
       }
     }
       
     else{
       return_val=MSQ_MAX_CAP;
     }
     
     return return_val;
   }

} //namespace


#endif // ConditionNumberQualityMetric_hpp


