// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file VolumeQualityMetric.hpp

Header file for the Mesquite::VolumeQualityMetric class

  \author Michael Brewer
  \date   2002-10-10
 */


#ifndef VolumeQualityMetric_hpp
#define VolumeQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityMetric.hpp"

namespace Mesquite
{
   class VolumeQualityMetric : public QualityMetric
  {
   public:
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~VolumeQualityMetric()
       {};

   protected:
 
    
   private:

    
  };


} //namespace


#endif // VolumeQualityMetric_hpp
