// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file CornerJacobianQualityMetric.hpp

Header file for the Mesquite::CornerJacobianQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef CornerJacobianQualityMetric_hpp
#define CornerJacobianQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "VolumeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"


namespace Mesquite
{
     /*! \class CornerJacobianQualityMetric
       \brief Computes the volume or area of the element, as appropriate.
       This metric uses the average of the corner Jacobian determinants
       for the approximation to the volume of hex.
     */
     /*!\todo MB:  CornerJacobianQualityMetric is currently not being
       evaluated using the corner jacobian method (instead it is using
       the area functions from MsqMeshEntity... needs to be modified
       to do so.*/
   class CornerJacobianQualityMetric : public VolumeQualityMetric
   {
  public:
 
       /*!Returns a pointer to a VolumeQualityMetric.  The metric
         does not use the sample point functionality or the
         compute_weighted_jacobian (except for possibly indirectly
         when evaluating the metric value for a hex).  It evaluates
         the signed area of surface elements and the signed volume
         of volume elements.  It does require a feasible region,
         and (in general) the metric needs to be minimized.
       */
     static VolumeQualityMetric* create_new(){
       
       VolumeQualityMetric* m = new CornerJacobianQualityMetric();
       return m;
     }
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~CornerJacobianQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,double &fval,
                           MsqError &err); 
          
  protected:

  private:
     
     CornerJacobianQualityMetric();
    
  };

} //namespace


#endif // CornerJacobianQualityMetric_hpp


