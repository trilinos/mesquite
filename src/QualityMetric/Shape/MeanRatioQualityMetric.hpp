// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file MeanRatioQualityMetric.hpp

Header file for the Mesquite::MeanRatioQualityMetric class

\author Michael Brewer
\author Thomas Leurent
\date   2002-06-19
 */


#ifndef MeanRatioQualityMetric_hpp
#define MeanRatioQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ShapeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
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
   private:
 
     MeanRatioQualityMetric() :
       ShapeQualityMetric()
     {
       MsqError err;
       set_metric_type(ELEMENT_BASED);
       set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
       set_negate_flag(1);
       set_gradient_type(ANALYTICAL_GRADIENT);
       set_hessian_type(ANALYTICAL_HESSIAN);
       avgMethod=QualityMetric::LINEAR;
       feasible=1;
       set_name("Mean Ratio");
     }

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
      bool evaluate_element(PatchData &pd, MsqMeshEntity *element, 
                            double &fval, MsqError &err); 

      bool compute_element_analytical_gradient(PatchData &pd,
                                               MsqMeshEntity *element,
                                               MsqVertex *free_vtces[], 
                                               Vector3D grad_vec[],
                                               int num_vtx, 
                                               double &metric_value,
                                               MsqError &err);

      bool compute_element_analytical_hessian(PatchData &pd,
                                              MsqMeshEntity *e,
                                              MsqVertex *v[], 
                                              Vector3D g[],
                                              Matrix3D h[],
                                              int nv, 
                                              double &m,
                                              MsqError &err);
   };
} //namespace


#endif // MeanRatioQualityMetric_hpp
