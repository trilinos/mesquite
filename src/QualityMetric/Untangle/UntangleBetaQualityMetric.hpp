// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file UntangleBetaQualityMetric.hpp

Header file for the Mesquite::UntangleBetaQualityMetric class

  \author Michael Brewer
  \date   2002-09-10
 */

#ifndef UNTANGLE_BETA_QUALITY_METRIC_HPP
#define UNTANGLE_BETA_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "UntangleQualityMetric.hpp"

namespace Mesquite
{
     /*! \class UntangleBetaQualityMetric
       Given a scalar value beta and local signed element volume alpha_i,
       define delta_i to be alpha_i minus beta.  The Untangle beta value
       is then defined as the sum over sample points of the absolute value
       of delta_i minus delta_i.
     */
   
   class UntangleBetaQualityMetric : public UntangleQualityMetric
   {
   public:
     
       /*! \fn UntangleQualityMetric* UntangleBetaQualityMetric::create_new()
         \brief The function create_new is used to create a untangle quality
         metric.  The constructor defaults to SUM AveragingMethod and
         ELEMENT_VERTICES evaluationMode.  The default beta value is
         .05.
       */
     static UntangleQualityMetric* create_new()
        {
          UntangleQualityMetric* m = new UntangleBetaQualityMetric();
          return m;
        }
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~UntangleBetaQualityMetric()
        {}
     
     double evaluate_element(PatchData &pd,
                             MsqMeshEntity* element,
                             MsqError &err);
     
   protected:
     
     
   private:
     double mBeta;
     UntangleBetaQualityMetric();
   };


} //namespace


#endif // UntangleBetaQualityMetric_hpp


