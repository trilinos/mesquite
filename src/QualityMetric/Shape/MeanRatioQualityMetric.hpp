// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file MeanRatioQualityMetric.hpp

Header file for the Mesquite::MeanRatioQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef MeanRatioQualityMetric_hpp
#define MeanRatioQualityMetric_hpp
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
   /*! \class MeanRatioQualityMetric
     \brief Computes the mean ratio of given element.
       
   */
   class MeanRatioQualityMetric : public ShapeQualityMetric
   {
   public:
 
      /*!Returns a pointer to a ShapeQualityMetric.  The metric
        does not use the sample point functionality or the
        compute_weighted_jacobian.  It evaluates the metric at
        the element vertices, and uses the isotropic ideal element.
        It does require a feasible region, and the metric needs
        to be minimized.
      */
      static ShapeQualityMetric* create_new() {
         ShapeQualityMetric* m = new MeanRatioQualityMetric();
         return m;
      }
     
      //! virtual destructor ensures use of polymorphism during destruction
      virtual ~MeanRatioQualityMetric() {
      }
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element, double &fval,
                             MsqError &err); 

  protected:     
     bool mean_ratio_2d(Vector3D temp_vec[],double &fval,MsqError &err);
     bool mean_ratio_3d(Vector3D temp_vec[],double &fval,MsqError &err);
  private:
     
     MeanRatioQualityMetric();
  };
     //BEGIN INLINE FUNCITONS
   inline bool MeanRatioQualityMetric::mean_ratio_2d(Vector3D temp_vec[],
                                                     double &fval,
                                                     MsqError &err)
   {
      // NOTE: Fixed to check for degenerate or inverted elements
      // NOTE: If coordinates are always first two, can make calculation faster
      double determinant = 1.0%(temp_vec[0]*temp_vec[1]);

      if (determinant <= MSQ_MIN) {
         return false;
      }

      fval = ((temp_vec[0].length_squared() +
               temp_vec[1].length_squared()) / (2.0*determinant));
      return true;
   }

   inline bool MeanRatioQualityMetric::mean_ratio_3d(Vector3D temp_vec[],
                                                     double &fval,
                                                       MsqError &err)
   {   
      // NOTE: Fixed to check for degenerate or inverted elements
      double determinant = temp_vec[0]%(temp_vec[1]*temp_vec[2]);

      if (determinant <= MSQ_MIN) {
         return false;
      }

      fval = ((temp_vec[0].length_squared() +
               temp_vec[1].length_squared() +
               temp_vec[2].length_squared()) / (3.0*pow(determinant,2.0/3.0)));
      return true;
   }

} //namespace


#endif // MeanRatioQualityMetric_hpp
