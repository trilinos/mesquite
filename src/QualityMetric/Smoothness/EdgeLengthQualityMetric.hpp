// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file EdgeLengthQualityMetric.hpp

Header file for the Mesquite::EdgeLengthQualityMetric class

  \author Michael Brewer
  \date   2002-06-13
 */


#ifndef EdgeLengthQualityMetric_hpp
#define EdgeLengthQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
     /*! \class EdgeLengthQualityMetric
       \brief Computes the lengths of the edges connected to given a vertex..
       
        EdgeLengthQualityMetric is a vertex based metric which computes
        the lengths of the edges connected to a given vertex and then
        averages those together, using the specified averaging method
        The metric uses SUM as the default averaging method.
     */
   class MsqMeshEntity;
   class MsqVertex;
   
   class EdgeLengthQualityMetric : public SmoothnessQualityMetric
  {
   public:
    
    EdgeLengthQualityMetric()
       {
         avgMethod=SUM;
         feasible=0;
         set_metric_type(QualityMetric::VERTEX_BASED);
         set_name("Edge Length Metric");
       }
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~EdgeLengthQualityMetric()
       {}

  
   protected:
   
    bool evaluate_vertex(PatchData &pd, MsqVertex *vert, double &fval,
                         MsqError &err);
  };


} //namespace


#endif // EdgeLengthQualityMetric_hpp


