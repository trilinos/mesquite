// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file UntangleQualityMetric.hpp

Header file for the Mesquite::UntangleQualityMetric class

  \author Michael Brewer
  \date   2002-10-10
 */


#ifndef UntangleQualityMetric_hpp
#define UntangleQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityMetric.hpp"

namespace Mesquite
{
   class UntangleQualityMetric : public QualityMetric
  {
   public:
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~UntangleQualityMetric()
       {};

   protected:
 
    
   private:

    
  };


} //namespace


#endif // UntangleQualityMetric_hpp
