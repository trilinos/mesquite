// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file PowerQualityMetric.hpp
\brief
Header file for the Mesquite::PowerQualityMetric class

  \author Michael Brewer
  \date   2002-09-05
 */


#ifndef PowerQualityMetric_hpp
#define PowerQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "CompositeQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
     /*! \class PowerQualityMetric
       \brief Raises a single quality metrics (qMetric1) to an arbitrary
       power (a double value, scaleAlpha) for two- and three-diminsional
       elements.  
     */
   class PowerQualityMetric : public CompositeQualityMetric
   {
  public:
       //! The function create_new is used to create a composite quality metric
       /*! Ensures that qm1 is not NULL.  If qm1 is only valid
         on a certain feasible, then the composite metric has the same
         constraint.  The composite metric also has the same negate flag
         as qm1.
       */
     static CompositeQualityMetric* create_new(QualityMetric* qm1,
                                               double pow_double, 
                                               MsqError &err){
      CompositeQualityMetric* m = new PowerQualityMetric(qm1, pow_double,
                                                         err);
      return m;
     }
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~PowerQualityMetric()
        {  }
     
     bool evaluate_element(PatchData& pd, MsqMeshEntity *element,double &value,
                             MsqError &err);
     bool evaluate_vertex(PatchData& pd, MsqVertex *vertex, double &value,
                            MsqError &err);
     
  protected:
     
  private:
     
     PowerQualityMetric(QualityMetric* qm1, double pow_double,MsqError &err);
     
   };
   

} //namespace


#endif // PowerQualityMetric_hpp







