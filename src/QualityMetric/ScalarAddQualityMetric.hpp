// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ScalarAddQualityMetric.hpp
\brief
Header file for the Mesquite::ScalarAddQualityMetric class

  \author Michael Brewer
  \date   April 14, 2003
 */


#ifndef ScalarAddQualityMetric_hpp
#define ScalarAddQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "CompositeQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
     /*! \class ScalarAddQualityMetric
       \brief Adds a number (a double) to the quality metric value.
     */
   class ScalarAddQualityMetric : public CompositeQualityMetric
   {
  public:
       /*! Ensures that qm1 is not NULL.  If qm1 is only valid
         on a certain feasible, then the composite metric has the same
         constraint.  The composite metric also has the same negate flag
         as qm1.
       */
     ScalarAddQualityMetric(QualityMetric* qm1, double scalar_double,
                            MsqError &err);
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~ScalarAddQualityMetric()
        {  }
     
     bool evaluate_element(PatchData& pd, MsqMeshEntity *element,double &value,
                             MsqError &err);
     bool evaluate_vertex(PatchData& pd, MsqVertex *vertex, double &value,
                            MsqError &err);
          
   };
   

} //namespace


#endif // ScalarAddQualityMetric_hpp







