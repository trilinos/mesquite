// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file MultiplyQualityMetric.hpp
\brief
Header file for the Mesquite::MultiplyQualityMetric class

  \author Michael Brewer
  \date   2002-09-05
 */


#ifndef MultiplyQualityMetric_hpp
#define MultiplyQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "CompositeQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
     /*! \class MultiplyQualityMetric
       \brief Combines two quality metrics (qMetric1 and qMetric2 defined
       in the parent class CompositeQualityMetric) by multiplication for two-
       and three-diminsional elements.  Note:  This function should not
       be used to combine a node-based metric with an element-based
       metric.  
     */
   class MultiplyQualityMetric : public CompositeQualityMetric
   {
  public:
       /*! Ensures that qm1 and qm2 are not NULL.  If either qm1 or qm2
         are valid only on a feasible region, then the composite
         metric's feasibility flag is set to one.  If qm1 and qm2 have
         different negateFlags, then a warning is printed, and the composite
         metric's negate flag is set to one.  Otherwise, the composite
         metric's negateFlag is set to qm1's negateFlag (and, thus, qm2's
         negateFlag).  
       */
     MultiplyQualityMetric(QualityMetric* qm1, QualityMetric* qm2,
                           MsqError &err);
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~MultiplyQualityMetric()
        {  }
     
     inline void set_multiplication_operands(QualityMetric* qm1,
                                             QualityMetric* qm2,
                                             MsqError &/*err*/){
       set_qmetric2(qm2);
       set_qmetric1(qm1);
     }

     bool evaluate_element(PatchData& pd, MsqMeshEntity *element,double &value,
                             MsqError &err);
     bool evaluate_vertex(PatchData& pd, MsqVertex *vertex, double &value,
                            MsqError &err);
     
  protected:
     
  private:
     
     
   };
   

} //namespace


#endif // MultiplyQualityMetric_hpp







