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
       
     */
   class ConditionNumberQualityMetric : public ShapeQualityMetric
   {
  public:
 
       /*!Returns a poitner to a ShapeQualityMetric.  The metric
         does not use the sample point functionality or the
         compute_weighted_jacobian.  It evaluates the metric at
         the element vertices, and uses the isotropic ideal element.
         It does require a feasible region, and the metric needs
         to be minimized.
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
          
  protected:     
     double condition_number_2d(Vector3D temp_vec[],MsqError &err);
     double condition_number_3d(Vector3D temp_vec[],MsqError &err);
  private:
     
     ConditionNumberQualityMetric();
    
  };
     //BEGIN INLINE FUNCITONS
   inline double ConditionNumberQualityMetric::condition_number_2d(Vector3D
                                                                   temp_vec[],
                                                                   MsqError
                                                                   &err)
   {
       //NOTE:  We take absolute value here
     double temp_val=fabs((temp_vec[0]*temp_vec[1]).length()*2.0);
     double return_val=MSQ_MAX_CAP;
     if(temp_val>MSQ_MIN){
       return_val=(temp_vec[0].length_squared()+temp_vec[1].length_squared())/
          temp_val;
     }
     return return_val;
   }

   inline double ConditionNumberQualityMetric::condition_number_3d(Vector3D
                                                                   temp_vec[],
                                                                   MsqError
                                                                   &err)
   {   
     double term1=temp_vec[0]%temp_vec[0]+
        temp_vec[1]%temp_vec[1]+
        temp_vec[2]%temp_vec[2];
       //norm squared of adjoint of J
     double term2=(temp_vec[0]*temp_vec[1])%
        (temp_vec[0]*temp_vec[1])+
        (temp_vec[1]*temp_vec[2])%
        (temp_vec[1]*temp_vec[2])+
        (temp_vec[2]*temp_vec[0])%
        (temp_vec[2]*temp_vec[0]);
       //det of J
     double temp_var=temp_vec[0]%(temp_vec[1]*temp_vec[2]);
     double return_val=sqrt(term1*term2);
     if(return_val>MSQ_MIN){
         //if not degenerate or inverted???
       return_val/=(3*temp_var);
     }
     else
       return_val=MSQ_MAX_CAP;
     return return_val;
   }
   
   

} //namespace


#endif // ConditionNumberQualityMetric_hpp


