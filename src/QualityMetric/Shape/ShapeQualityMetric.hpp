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
      /*! \class ShapeQualityMetric
       \brief Parent class for the Shape Quality Metrics.
       
     */
   public:
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeQualityMetric()
       {};

   protected:
 
      //void compute_scalar_weights(int num_scalar_weights, double scalar_weights[], MsqError &err);

      //! Given the 2-d jacobian matrix, compute the condition number, fval 
    bool condition_number_2d(Vector3D temp_vec[],size_t v_ind, PatchData &pd,
                              double &fval, MsqError &err);
      //! Given the 3-d jacobian matrix, compute the condition number, fval 
    bool condition_number_3d(Vector3D temp_vec[], PatchData &pd, double &fval, MsqError &err);
    
   private:

    
  };

   /*
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
       Vector3D unit_surf_norm;
       pd.get_domain_normal_at_vertex(v_ind,true,unit_surf_norm,err);MSQ_CHKERR(err);
       //if invalid 
       if(unit_surf_norm%cross_vec < 0.0){
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
   */

#undef __FUNC__
#define __FUNC__ "ShapeQualityMetric::condition_number_2d"
  inline bool ShapeQualityMetric::condition_number_2d(Vector3D temp_vec[],
                                                       size_t v_ind,
                                                       PatchData &pd,
                                                       double &fval,
                                                       MsqError &err)
   {   
       //norm squared of J
     double term1=temp_vec[0]%temp_vec[0]+
        temp_vec[1]%temp_vec[1];

     Vector3D unit_surf_norm;
     if ( pd.domain_set() ) {
       pd.get_domain_normal_at_vertex(v_ind,true,unit_surf_norm,err);MSQ_CHKERR(err);
     }
     else {
        return false;
     }

     // det J
     double temp_var=unit_surf_norm%(temp_vec[0]*temp_vec[1]);

     double h;
     double delta=pd.get_barrier_delta_2d(err); 
     // Note: technically, we want delta=eta*tau-max
     //       whereas the function above gives delta=eta*alpha-max
     //      
     //       Because the only requirement on eta is eta << 1,
     //       and because tau-max = alpha-max/0.707 we can
     //       ignore the discrepancy

     if (delta==0) { 
        h=temp_var;
     }
     else {
        h = vertex_barrier_function(temp_var,delta);
     }
     // Note: when delta=0, the vertex_barrier_function
     //       formally gives h=temp_var as well.
     //       We just do it this way to avoid any 
     //       roundoff issues.
     // Also: when delta=0, this metric is identical
     //       to the original condition number with
     //       the barrier at temp_var=0

     if (h<MSQ_DBL_MIN) {
       err.set_msg("Barrier function is zero due to excessively large negative area compared to delta. /n Try to untangle mesh another way. ");
       return false;
     }

     fval=term1/(2*h);

     if (fval>MSQ_MAX_CAP) {
        return false; 
     }
     return true;

   }
   

   //} //namespace

#undef __FUNC__
#define __FUNC__ "ShapeQualityMetric::condition_number_3d"
  inline bool ShapeQualityMetric::condition_number_3d(Vector3D temp_vec[],
                                                       PatchData &pd,
                                                       double &fval,
                                                       MsqError &err)
   {   
       //norm squared of J
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

     double h;
     double delta=pd.get_barrier_delta_3d(err); 
     // Note: technically, we want delta=eta*tau-max
     //       whereas the function above gives delta=eta*alpha-max
     //      
     //       Because the only requirement on eta is eta << 1,
     //       and because tau-max = alpha-max/0.707 we can
     //       ignore the discrepancy

     if (delta==0) { 
        h=temp_var;
     }
     else {
        h = vertex_barrier_function(temp_var,delta);
     }
     // Note: when delta=0, the vertex_barrier_function
     //       formally gives h=temp_var as well.
     //       We just do it this way to avoid any 
     //       roundoff issues.
     // Also: when delta=0, this metric is identical
     //       to the original condition number with
     //       the barrier at temp_var=0

     if (h<MSQ_DBL_MIN) {
       err.set_msg("Barrier function is zero due to excessively large negative area compared to delta. /n Try to untangle mesh another way. ");
       return false;
     }

     fval=sqrt(term1*term2)/(3*h);

     if (fval>MSQ_MAX_CAP) {
        return false; 
     }
     return true;

   }
   

} //namespace


#endif // ShapeQualityMetric_hpp
