// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file InverseMeanRatioQualityMetric.hpp

Header file for the Mesquite::InverseMeanRatioQualityMetric class

  \author Michael Brewer
  \date   2002-11-11
 */


#ifndef InverseMeanRatioQualityMetric_hpp
#define InverseMeanRatioQualityMetric_hpp
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
     /*! \class InverseMeanRatioQualityMetric
       \brief Computes the inverse mean ratio quality metric
       of given element.
       
     */
   class InverseMeanRatioQualityMetric : public ShapeQualityMetric
   {
  public:
 
       /*!Returns a pointer to a ShapeQualityMetric.  The metric
         does not use the sample point functionality or the
         compute_weighted_jacobian.  It evaluates the metric at
         the element vertices, and uses the isotropic ideal element.
         It does require a feasible region, and the metric needs
         to be maximized.
       */
     static ShapeQualityMetric* create_new(){
       
       ShapeQualityMetric* m = new InverseMeanRatioQualityMetric();
       return m;
     }
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~InverseMeanRatioQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     double evaluate_element(PatchData &pd, MsqMeshEntity *element,
                             MsqError &err); 
          
  protected:     
     double inverse_mean_ratio_2d(Vector3D temp_vec[],MsqError &err);
     double inverse_mean_ratio_3d(Vector3D temp_vec[],MsqError &err);
  private:
     
     InverseMeanRatioQualityMetric();
    
  };
     //BEGIN INLINE FUNCITONS
   inline double InverseMeanRatioQualityMetric::inverse_mean_ratio_2d(Vector3D
                                                                   temp_vec[],
                                                                   MsqError
                                                                   &err)
   {
       //NOTE:  We take absolute value here
     double area_val=fabs((temp_vec[0]*temp_vec[1]).length()*2.0);
     double return_val=(temp_vec[0].length_squared()+
                        temp_vec[1].length_squared());
     if(return_val>MSQ_MIN){
       return_val=area_val/return_val;
     }
     else{
       return_val=0.0;
     }
     return return_val;
   }

   inline double InverseMeanRatioQualityMetric::inverse_mean_ratio_3d(Vector3D
                                                                   temp_vec[],
                                                                   MsqError
                                                                   &err)
   {   
     double term1=temp_vec[0]%temp_vec[0]+
        temp_vec[1]%temp_vec[1]+
        temp_vec[2]%temp_vec[2];
       //det of J
     double temp_vol_sqr=temp_vec[0]%(temp_vec[1]*temp_vec[2]);
     temp_vol_sqr*=temp_vol_sqr;
     double return_val=0.0;
     if(term1>MSQ_MIN){
         //if not degenerate or inverted???
       return_val=3*(pow(temp_vol_sqr,(1.0/3.0)))/term1;
     }
     return return_val;
   }
   
   

} //namespace


#endif // InverseMeanRatioQualityMetric_hpp


