/*! \file VertexConditionNumberQualityMetric.hpp

Header file for the Mesquite::VertexConditionNumberQualityMetric class

  \author Michael Brewer
  \date   April 14, 2003
 */


#ifndef VertexConditionNumberQualityMetric_hpp
#define VertexConditionNumberQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ShapeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"


namespace Mesquite
{
     /*! \class VertexConditionNumberQualityMetric
       \brief Computes the condition numbers of the corner's of elements
       connected to the given vertex and then averages those values.

       The metric does not use the sample point functionality or the
       compute_weighted_jacobian.  It uses the isotropic ideal
       element.  This metric does require a feasible region, and
       the metric needs to be minimized.
     */
   class VertexConditionNumberQualityMetric : public ShapeQualityMetric
   {
  public:
     VertexConditionNumberQualityMetric();
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~VertexConditionNumberQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_vertex(PatchData &pd, MsqVertex *vert, double &fval,
                           MsqError &err); 
          
  };
    
   

} //namespace

#endif // VertexConditionNumberQualityMetric_hpp

