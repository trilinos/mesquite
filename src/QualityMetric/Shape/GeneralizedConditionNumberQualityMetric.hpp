// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file GeneralizedConditionNumberQualityMetric.hpp

Header file for the Mesquite::GeneralizedConditionNumberQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef GeneralizedConditionNumberQualityMetric_hpp
#define GeneralizedConditionNumberQualityMetric_hpp
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
     /*! \class GeneralizedConditionNumberQualityMetric
       \brief Computes the condition number of given element.
       The``condition number" is scaled between one and infinity,
       with an ideal element having condition number one.  
     */
   class GeneralizedConditionNumberQualityMetric : public ShapeQualityMetric
   {
  public:
 
       /*!Returns a pointer to a ShapeQualityMetric.  Defaults the metric
         to the LINEAR averaging method and ELEMENT_VERTICES sample points.
         The default metric name is "Condition Number".  It does require a
         feasible region, and the metric needs to be maximized.
       */
     static ShapeQualityMetric* create_new(){
       
       ShapeQualityMetric* m = new GeneralizedConditionNumberQualityMetric();
       return m;
     }
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~GeneralizedConditionNumberQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element, double &fval,
                           MsqError &err);
     
     //! Evaluate the "condition number" for a vertex
     bool evaluate_vertex(PatchData &pd, MsqVertex *vertex, double &fval,
                          MsqError &err);
     
  protected:

     bool compute_condition_number(PatchData &pd, MsqMeshEntity *elem,
                                   Vector3D* jacobian_vectors,
                                   int num_jacobian_vectors,
                                   double &fval,
                                   MsqError &err);
     
     
  private:
     
     GeneralizedConditionNumberQualityMetric();
    
  };
//BEGIN INLINE FUNCTIONS

   inline bool GeneralizedConditionNumberQualityMetric::compute_condition_number(
      PatchData &pd, MsqMeshEntity *element, Vector3D* jacobian_vectors,
      int num_jacobian_vectors, double &fval, MsqError &err)
   {
       //PRINT_INFO("INSIDE CONDITION NUMBER COMPUTE_CON\n");
     double temp_var=0;
     if(num_jacobian_vectors==2){
       size_t vert=element->get_vertex_index(0);
       Vector3D cross_vec=jacobian_vectors[0]*jacobian_vectors[1];
       if ( pd.domain_set() ) {
         Vector3D norm_vec;
         pd.get_domain_normal_at_vertex(vert,norm_vec,err);MSQ_CHKERR(err);
         if(cross_vec%norm_vec<0.0){
           return false;
         }
       }

       temp_var=fabs((cross_vec).length());
       fval=jacobian_vectors[0].length_squared();
       fval+=jacobian_vectors[1].length_squared();
       if(temp_var>=MSQ_MIN){ //if not degenerate
         fval/=(2.0*temp_var);
       }
       else{
         fval=MSQ_MAX_CAP;
       }
       return true;
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
       if(temp_var<=0.0){
         return false;
       }
       
       fval=sqrt(term1*term2);
       if(fval>MSQ_MIN){
           //if not degenerate or inverted???
         fval/=(3*temp_var);
       }
       else{
         fval=MSQ_MAX_CAP;
       }
     }
       
     else{
       fval=MSQ_MAX_CAP;
     }
     
     return true;
   }

} //namespace


#endif // GeneralizedConditionNumberQualityMetric_hpp


