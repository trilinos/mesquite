// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ASMQualityMetric.hpp

Header file for the Mesquite::ASMQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef ASMQualityMetric_hpp
#define ASMQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"


namespace Mesquite
{
     /*! \class ASMQualityMetric
       \brief Computes the ASM (Area Smoothness Metric) of the given element.
       
        The metric does not use the sample point functionality or the
        compute_weighted_jacobian.  It computes the unsigned area
        or volume (a_0)of the element and of neighboring elements (a_i).
        In 2-d, elements are considered neighbors if they share an
        edge.  In 3-d, elements are considered neighbors if they
        share a face.  The metric is then taken to be the average
        (default is MAXIMUM) of abs(a_i-a_0)/(a_i+a_0) , and the metric
        needs to be minimized.  The ideal metic value is zero which occurs
        when an element and its neighbors have the same unsigned area
        (or volume in three-dimensions).
     */
   class ASMQualityMetric : public SmoothnessQualityMetric
   {
  public:
 
       /*!Returns a pointer to a SmoothnessQualityMetric. 
         
       */
     static SmoothnessQualityMetric* create_new(){
       
       SmoothnessQualityMetric* m = new ASMQualityMetric();
       return m;
     }
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~ASMQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,double &fval,
                           MsqError &err); 
          
  protected:
  private:
     
     ASMQualityMetric();
    
  };

   
   

} //namespace


#endif // ASMQualityMetric_hpp


