// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file EdgeLengthRangeQualityMetric.hpp

Header file for the Mesquite::EdgeLengthRangeQualityMetric class

  \author Michael Brewer
  \date   2002-06-13
 */


#ifndef EdgeLengthRangeQualityMetric_hpp
#define EdgeLengthRangeQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
     /*! \class EdgeLengthRangeQualityMetric
       \brief Computes the edge length range metric for a given vertex.
       
        EdgeLengthRangeQualityMetric is a vertex based metric which computes
        the lengths of the edges connected to a given vertex and then
        uses those values to form a metric.  The metric is created using
        two doubles A and B.  The value of the metric is zero (ideal) if
        the edge lengths are in the range [A,B].  Otherwise, the
        metric value is some positive number.  For a given vertex,
        v_i, with connected edges of lengths l_j for j=1...k, the metric
        value is the average (where the default average type is SUM) of
        u_j = ( | l_j - A | - (l_j - A) )^2 + ( | B - l_j | - (B - l_j) )^2.
     */
   class MsqMeshEntity;
   class MsqVertex;
   
   class EdgeLengthRangeQualityMetric : public SmoothnessQualityMetric
  {
   public:
    
      //This is the form of the constructor that should beused.
    EdgeLengthRangeQualityMetric(double low_a, double high_a, MsqError &err)
       {
         if(low_a>high_a){
           err.set_msg("Edge Length Range values given in descending order.");
         }
         lowVal=low_a;
         highVal=high_a;
         avgMethod=SUM;
         feasible=0;
         set_metric_type(QualityMetric::VERTEX_BASED);
         set_name("Edge Length Range Metric");
       }

      // virtual destructor ensures use of polymorphism during destruction
    virtual ~EdgeLengthRangeQualityMetric()
       {}

  
   protected:
   
    bool evaluate_vertex(PatchData &pd, MsqVertex *vert, double &fval,
                         MsqError &err);

   private:
    
      //Generally, this constructor should not be used.
    EdgeLengthRangeQualityMetric()
       {
         lowVal=0;
         highVal=0;
         avgMethod=SUM;
         feasible=0;
         set_metric_type(QualityMetric::VERTEX_BASED);
         set_name("Edge Length Range Metric Default Constructor");
       }
    
    double highVal;
    double lowVal;
    
  };


} //namespace


#endif // EdgeLengthRangeQualityMetric_hpp


