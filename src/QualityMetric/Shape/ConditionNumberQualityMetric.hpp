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

       The metric does not use the sample point functionality or the
       compute_weighted_jacobian.  It evaluates the metric at
       the element vertices, and uses the isotropic ideal element.
       It does require a feasible region, and the metric needs
       to be minimized.
     */
   class ConditionNumberQualityMetric : public ShapeQualityMetric
   {
  public:
     ConditionNumberQualityMetric();
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~ConditionNumberQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,double &fval,
                           MsqError &err); 
          
  protected:
     
  private:
    
  };
    
   

} //namespace


#endif // ConditionNumberQualityMetric_hpp


