// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file SmoothnessQualityMetric.hpp

Header file for the Mesquite::SmoothnessQualityMetric class

  \author Thomas Leurent
  \date   2002-09-01
 */


#ifndef SmoothnessQualityMetric_hpp
#define SmoothnessQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityMetric.hpp"

namespace Mesquite
{
   class MsqMeshEntity;
   class SmoothnessQualityMetric : public QualityMetric
  {
   public:
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~SmoothnessQualityMetric()
       {};

   protected:
 
    void compute_scalar_weights(int num_scalar_weights, double scalar_weights[], MsqError &err);

   private:

    
  };


} //namespace


#endif // SmoothnessQualityMetric_hpp
