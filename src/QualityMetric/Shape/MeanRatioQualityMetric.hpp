// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file MeanRatioQualityMetric.hpp

Header file for the Mesquite::MeanRatioQualityMetric class

  \author Michael Brewer
  \date   2002-09-05
 */

#ifndef MEAN_RATIO_QUALITY_METRIC_HPP
#define MEAN_RATIO_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "ShapeQualityMetric.hpp"

namespace Mesquite
{
     /*! \class MeanRatioQualityMetric
       \brief MeanRatioQualityMetric evaluates the mean ratio
       for a given elemnt.
     */
   
   class MeanRatioQualityMetric : public ShapeQualityMetric
   {
   public:
     
       /*! \fn ShapeQualityMetric* MeanRatioQualityMetric::create_new()
         \brief The function create_new is used to create a shape quality
         metric.  The constructor defaults to LINEAR AveragingMethod and
         ELEMENT_VERTICES evaluationMode.
       */
     static ShapeQualityMetric* create_new()
        {
          ShapeQualityMetric* m = new MeanRatioQualityMetric();
          return m;
        }
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~MeanRatioQualityMetric()
        {}
     
     double evaluate_element(PatchData &pd,
                             MsqMeshEntity* element,
                             MsqError &err);
     
   protected:
     
     
   private:
     
     MeanRatioQualityMetric();
   };


} //namespace


#endif // MeanRatioQualityMetric_hpp


