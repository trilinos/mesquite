/*!
  \file   LInfTemplate.cpp
  \brief  

  This Objective Function is evaluated using an L-infinity norm.
  total=max (abs(x))
  \author Michael Brewer
  \date   2002-07-3
*/
#include <math.h>
#include "LInfTemplate.hpp"
using  namespace Mesquite;  

#undef __FUNC__
#define __FUNC__ "LInfTemplate::LInfTemplate"

LInfTemplate::LInfTemplate(QualityMetric *qualitymetric){
   set_quality_metric(qualitymetric);
   set_feasible(qualitymetric->get_feasible_constraint());
   set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
   set_negate_flag(qualitymetric->get_negate_flag());
}

#undef __FUNC__
#define __FUNC__ "LInfTemplate::~LInfTemplate"

//Michael:  need to clean up here
LInfTemplate::~LInfTemplate(){

}

#undef __FUNC__
#define __FUNC__ "LInfTemplate::concrete_evaluate"
bool LInfTemplate::concrete_evaluate(PatchData &patch, double &fval,
                                     MsqError &err){
    //Total value of objective function
  double temp_value=0;
  fval=0.0;
  bool obj_bool=true;
    //For elements in Patch
  int index;
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL)
    currentQM=get_quality_metric_list().front();
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED){
    int num_elements=patch.num_elements();
    MsqMeshEntity* elems=patch.get_element_array(err);
      //Michael:  this may not do what we want
      //Set currentQM to be the first quality metric* in the list 
    for (index=0; index<num_elements;index++){
        //evaluate metric for this elem
      obj_bool=currentQM->evaluate_element(patch, &elems[index], temp_value,
                                           err);
      MSQ_CHKERR(err);
        //if invalid patch
      if(!obj_bool){
        return false;
      }
      temp_value = fabs(temp_value);
      MSQ_CHKERR(err);
      if(temp_value>fval)
        fval=temp_value;
    }//end loop over elements
  }//end if not VERTEX
  else {//VERTEX
    int num_vertices=patch.num_vertices();
    MsqVertex* vertices=patch.get_vertex_array(err);
      //Michael:  this may not do what we want
      //Set currentQM to be the first quality metric* in the list
 
    for (index=0; index<num_vertices;index++){
        //evaluate metric for this vertex
      obj_bool=currentQM->evaluate_vertex(patch,&vertices[index],
                                          temp_value, err);
      MSQ_CHKERR(err);
        //if invalid patch
      if(!obj_bool){
        return false;
      }
      temp_value = fabs(temp_value);
      if(temp_value>fval)
        fval=temp_value;
    }//end loop over vertices
  }//end elseVERTEX
  return true;
}
	
	
