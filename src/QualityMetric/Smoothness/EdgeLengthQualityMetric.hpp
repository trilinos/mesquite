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
   class MsqMeshEntity;
   class MsqVertex;
   
   class EdgeLengthQualityMetric : public SmoothnessQualityMetric
  {
   public:
      // The function create_new is used to create a shape quality metric
      static SmoothnessQualityMetric* create_new(){

        EdgeLengthQualityMetric* m = new EdgeLengthQualityMetric();
        return m;
      }
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~EdgeLengthQualityMetric()
       {}

  
   protected:
   
    bool evaluate_vertex(PatchData &pd, MsqVertex *vert, double &fval,
                         MsqError &err);

   private:

    EdgeLengthQualityMetric()
       {
         avgMethod=SUM;
         feasible=0;
         set_metric_type(QualityMetric::VERTEX_BASED);
         set_name("Edge Length Metric");
       }
    
  };


} //namespace


#endif // EdgeLengthQualityMetric_hpp


