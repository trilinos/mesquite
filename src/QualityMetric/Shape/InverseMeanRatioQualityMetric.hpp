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
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,
                           double &fval, MsqError &err); 
          
  protected:     
     bool inverse_mean_ratio_2d(Vector3D temp_vec[],size_t v_ind,
                                PatchData &pd, double &fval,
                                MsqError &err);
     bool inverse_mean_ratio_3d(Vector3D temp_vec[],double &fval,
                                MsqError &err);
  private:
     
     InverseMeanRatioQualityMetric();
    
  };
     //BEGIN INLINE FUNCITONS
   inline bool InverseMeanRatioQualityMetric::inverse_mean_ratio_2d(Vector3D
                                                                    temp_vec[],
                                                                    size_t v_ind,
                                                                    PatchData &pd,
                                                                    double
                                                                    &fval,
                                                                    MsqError
                                                                    &err)
   {
     Vector3D cross_vec=(temp_vec[0]*temp_vec[1]);
     double area_val=fabs((cross_vec).length()*2.0);
       //NOTE:: the equal below is ONLY to get the code to work
       //when no normal is available.
     Vector3D surf_norm=cross_vec;
     pd.get_surface_normal(v_ind,surf_norm,err);MSQ_CHKERR(err);
       //if invalid
     if(surf_norm%cross_vec < 0.0){
       return false;
     }
     fval=(temp_vec[0].length_squared()+temp_vec[1].length_squared());
     if(fval>MSQ_MIN){
       fval=area_val/fval;
     }
     else{
       fval=0.0;
     }
     return true;
   }

   inline bool InverseMeanRatioQualityMetric::inverse_mean_ratio_3d(Vector3D
                                                                    temp_vec[],
                                                                    double
                                                                    &fval,
                                                                    MsqError
                                                                    &err)
   {   
     double term1=temp_vec[0]%temp_vec[0]+
        temp_vec[1]%temp_vec[1]+
        temp_vec[2]%temp_vec[2];
       //det of J
     double temp_vol_sqr=temp_vec[0]%(temp_vec[1]*temp_vec[2]);
       //if inverted
     if(temp_vol_sqr<=0.0)
       return false;
     
     temp_vol_sqr*=temp_vol_sqr;
     fval=0.0;
     if(term1>MSQ_MIN){
         //if not degenerate (or inverted)
       fval=3*(pow(temp_vol_sqr,(1.0/3.0)))/term1;
     }
     return true;
   }
   
   

} //namespace


#endif // InverseMeanRatioQualityMetric_hpp


