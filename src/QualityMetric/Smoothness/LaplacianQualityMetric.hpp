// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file LaplacianQualityMetric.hpp

Header file for the Mesquite::LaplacianQualityMetric class

  \author Michael Brewer
  \date   2002-06-13
 */


#ifndef LaplacianQualityMetric_hpp
#define LaplacianQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
   class MsqMeshEntity;
   class MsqVertex;
   
   class LaplacianQualityMetric : public SmoothnessQualityMetric
  {
   public:
      // The function create_new is used to create a shape quality metric
      static SmoothnessQualityMetric* create_new(){

        LaplacianQualityMetric* m = new LaplacianQualityMetric();
        return m;
      }
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~LaplacianQualityMetric()
       {}

  
   protected:
   
    double evaluate_node(MsqVertex *node, MsqError &err);

   private:

    LaplacianQualityMetric()
       {
         avgMethod=SUM;
         feasible=0;
         evalMode=QualityMetric::VERTEX;
         set_name("Laplacian Metric");
       }
    
  };


} //namespace


#endif // LaplacianQualityMetric_hpp


