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


namespace Mesquite
{
     /*! \class ConditionNumberQualityMetric
       \brief Computes the condition number of given element.
       
     */
   class ConditionNumberQualityMetric : public ShapeQualityMetric
   {
  public:
 
       /*!Returns a pointer to a ShapeQualityMetric.  The metric
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
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,double &fval,
                           MsqError &err); 
          
  protected:
     bool condition_number_2d(Vector3D temp_vec[],double &fval, MsqError &err);
     bool condition_number_3d(Vector3D temp_vec[],double &fval, MsqError &err);
  private:
     
     ConditionNumberQualityMetric();
    
  };
     //BEGIN INLINE FUNCITONS
   inline bool ConditionNumberQualityMetric::condition_number_2d(Vector3D
                                                                   temp_vec[],
                                                                   double
                                                                   &fval,
                                                                   MsqError
                                                                   &err)
   {
       //NOTE:  We take absolute value here
     double temp_val=fabs((temp_vec[0]*temp_vec[1]).length()*2.0);
     fval=MSQ_MAX_CAP;
     if(temp_val>MSQ_MIN){
       fval=(temp_vec[0].length_squared()+temp_vec[1].length_squared())/
          temp_val;
     }
       //returning true always until surf_normal is avail.
     return true;
   }

   inline bool ConditionNumberQualityMetric::condition_number_3d(Vector3D
                                                                 temp_vec[],
                                                                 double &fval,
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
       //if invalid
     if(temp_var<=0.0){
       fval=MSQ_MAX_CAP;
       return false;
     }
     
     fval=sqrt(term1*term2);
     if(fval>MSQ_MIN){
         //if not degenerate
       fval/=(3*temp_var);
     }
     else
       fval=MSQ_MAX_CAP;
       //return true
     return true;
   }
   
   

} //namespace


#endif // ConditionNumberQualityMetric_hpp


