/*!
  \file   MaxTemplate.cpp
  \brief  

  This Objective Function is the maximum of the quality metrics
  total = max(x)
  \author Lori Freitag
  \date   2002-07-18
*/
#include <math.h>
#include "MaxTemplate.hpp"
using  namespace Mesquite;  

#undef __FUNC__
#define __FUNC__ "MaxTemplate::MaxTemplate"

MaxTemplate::MaxTemplate(QualityMetric *qualitymetric){
   set_quality_metric(qualitymetric);
   set_negate_flag(qualitymetric->get_negate_flag());
}

#undef __FUNC__
#define __FUNC__ "MaxTemplate::~MaxTemplate"

//Lori:  need to clean up here
MaxTemplate::~MaxTemplate(){

}

#undef __FUNC__
#define __FUNC__ "MaxTemplate::concrete_evaluate"

bool MaxTemplate::concrete_evaluate(PatchData &patch, double &fval,
                                    MsqError &err){

  //Total value of objective function
  fval = 0.0;
  double temp_value=0;
  bool obj_bool = true;
  //For elements in Patch
  int index;
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED){
    
    int num_elements=patch.num_elements();
    MsqMeshEntity* elems=patch.get_element_array(err);
 
    for (index=0; index<num_elements; index++){

      //evaluate metric for this elem
      obj_bool = currentQM->evaluate_element(patch, &elems[index],
                                             temp_value, err);
      MSQ_CHKERR(err);
        //if invalid patch
      if(! obj_bool ){
        fval = 0.0;
        return false;
      }

      if(temp_value>fval )
        fval=temp_value;

    }//end loop over elements
  }//end if not VERTEX

  else if (currentQM->get_metric_type()==QualityMetric::VERTEX_BASED) {

    int num_vertices=patch.num_vertices();
    MsqVertex* vertices=patch.get_vertex_array(err);
 
    for (index=0; index<num_vertices;index++){

      //evaluate metric for this vertex
      obj_bool=currentQM->evaluate_vertex(patch, &vertices[index],
                                          temp_value, err);
      MSQ_CHKERR(err);
      //if invalid patch
      if(! obj_bool ){
        fval = 0.0;
        return false;
      }
      if(temp_value>fval)
        fval=temp_value;

    }//end loop over vertices
  }//end else VERTEX
  else {
    err.set_msg("Make sure MetricType is initialised in concrete QualityMetric constructor.");
  }
  
  return true;
}
	
	
