// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file AspectRatioGammaQualityMetric.hpp
  \brief
  Header file for the Mesquite::AspectRatioGammaQualityMetric class

  \author Michael Brewer
  \date   2002-05-16
 */


#ifndef AspectRatioGammaQualityMetric_hpp
#define AspectRatioGammaQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ShapeQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
     /*! \class AspectRatioGammaQualityMetric
       \brief Object for computing the aspect ratio gamma of
       simplicial elements.
     */
   class AspectRatioGammaQualityMetric : public ShapeQualityMetric
   {
   public:
       /*! The function create_new is used to create a shape quality metric.
        */
     static ShapeQualityMetric* create_new()
        {
          ShapeQualityMetric* m = new AspectRatioGammaQualityMetric();
          return m;
        }
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~AspectRatioGammaQualityMetric()
        {}
     
   protected:
     
     
   private:
     
     AspectRatioGammaQualityMetric()
        {
          feasible=0;
          evalMode=QualityMetric::ELEMENT_VERTICES;
          set_name("Aspect Ratio Gamma");
        }
       //!Returns the aspect ratio gamma of element.  If element
       //!is not a tetrahedron or triangle, sets an error.
     double evaluate_element(PatchData& pd,
                             MsqMeshEntity* element,
                             MsqError &err);
   };
   
   
} //namespace


#endif // AspectRatioGammaQualityMetric_hpp


