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

namespace Mesquite
{
   class ShapeQualityMetric : public QualityMetric
  {
   public:
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeQualityMetric()
       {};

   protected:
 
    void compute_scalar_weights(int num_scalar_weights, double scalar_weights[], MsqError &err);

   private:

    
  };


} //namespace


#endif // ShapeQualityMetric_hpp
