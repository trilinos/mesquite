// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ShapeQualityMetric.hpp

Header file for the Mesquite::ShapeQualityMetric class

  \author Thomas Leurent
  \date   2002-09-01
 */


#ifndef ShapeQualityMetric_hpp
#define ShapeQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityMetric.hpp"
#include "PatchData.hpp"

namespace Mesquite
{
   class ShapeQualityMetric : public QualityMetric
  {
   public:
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeQualityMetric()
       {};

   protected:
 
      //void compute_scalar_weights(int num_scalar_weights, double scalar_weights[], MsqError &err);

      //given the 2-d jacobian matrix, compute the condition number, fval 
    bool condition_number_2d(Vector3D temp_vec[],size_t v_ind, PatchData &pd,
                              double &fval, MsqError &err);
      //given the 3-d jacobian matrix, compute the condition number, fval 
    bool condition_number_3d(Vector3D temp_vec[],double &fval, MsqError &err);
    
   private:

    
  };


  //BEGIN INLINE FUNCITONS
#undef __FUNC__
#define __FUNC__ "ShapeQualityMetric::condition_number_2d"
   inline bool ShapeQualityMetric::condition_number_2d(Vector3D temp_vec[],
                                                       size_t v_ind,
                                                       PatchData &pd,
                                                       double &fval,
                                                       MsqError &err)
   {       
     Vector3D cross_vec=(temp_vec[0]*temp_vec[1]);
       //If the domain is not set, we assume all elements are valid.
       //Otherwise, we ensure the surface normal and the cross
       //vector have generally the same direction.
     if ( pd.domain_set() ) {
       Vector3D surf_norm;
       pd.get_domain_normal_at_vertex(v_ind,surf_norm,err);MSQ_CHKERR(err);
       //if invalid 
       if(surf_norm%cross_vec < 0.0){
         return false;
       }
     }
     
     double temp_val=cross_vec.length()*2.0;
     fval=MSQ_MAX_CAP;
     if(temp_val>MSQ_MIN){
       fval=(temp_vec[0].length_squared()+temp_vec[1].length_squared())/
          temp_val;
     }
       //returning true always until surf_normal is avail.
     return true;
   }


#undef __FUNC__
#define __FUNC__ "ShapeQualityMetric::condition_number_3d"
  inline bool ShapeQualityMetric::condition_number_3d(Vector3D temp_vec[],
                                                       double &fval,
                                                       MsqError &/*err*/)
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


#endif // ShapeQualityMetric_hpp
