/*!
  \file   LPTemplate.cpp
  \brief  

  This Objective Function is evaluated using an L P norm.
  total=(sum (x_i)^pVal)^(1/pVal)
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "LPTemplate.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
using  namespace Mesquite;  

#undef __FUNC__
#define __FUNC__ "LPTemplate::LPTemplate"

LPTemplate::LPTemplate(QualityMetric *qualitymetric, int Pinput, MsqError &err){
  set_quality_metric(qualitymetric);
  pVal=Pinput;
  if(pVal<2){
    err.set_msg("P_VALUE must be greater than 1.");
  }
  set_feasible(qualitymetric->get_feasible_constraint());
  set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
  set_negate_flag(qualitymetric->get_negate_flag());
}

#undef __FUNC__
#define __FUNC__ "LPTemplate::~LPTemplate"

//Michael:  need to clean up here
LPTemplate::~LPTemplate(){

}

#undef __FUNC__
#define __FUNC__ "LPTemplate::concrete_evaluate"
double LPTemplate::concrete_evaluate(PatchData &patch, MsqError &err){
    //Total value of objective function
//  double total_value=0;
//  double temp_value=0;
  int index=0;
//  double accum=0;
  MsqMeshEntity* elems=patch.get_element_array(err);
  
    //double check for pVal=0;
  if(pVal==0){
    err.set_msg("pVal equal zero not allowed.  L_0 is not a valid norm.");
    return 0;
  }
  
    //Michael:  this may not do what we want
    //Set currentQM to be the first quality metric* in the list
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL)
    currentQM=get_quality_metric_list().front();
//  MsqMeshEntity* current_ent;
  int num_elements=patch.num_elements();
  int num_vertices=patch.num_vertices();
  int total_num=0;
  if(currentQM->get_evaluation_mode()!=QualityMetric::VERTEX)   
    total_num=num_elements;
  else
    total_num=num_vertices;
  double *metric_values= new double[total_num];
  if(currentQM->get_evaluation_mode()!=QualityMetric::VERTEX)
  {
//    MsqMeshEntity* current_ent;
    for (index=0; index<num_elements;index++)
    {
      metric_values[index]=fabs(currentQM->evaluate_element(patch,
                                                            (&elems[index]),
                                                            err));
      MSQ_DEBUG_ACTION(3,{std::cout<< "      o  Quality metric value for element "
                          << index << "\t: " << metric_values[index] << "\n";});
    }
  }
  else
  {
    MsqVertex* vertices=patch.get_vertex_array(err);
    for (index=0; index<num_vertices;index++)
    {
        //evaluate metric for this vertex
      metric_values[index]=fabs(currentQM->evaluate_vertex(patch,
                                                           (&vertices[index]),
                                                           err));
    }
  }
  double obj_val=compute_function(metric_values, total_num, err);
  delete metric_values;
  return obj_val;
}

#undef __FUNC__
#define __FUNC__ "LPTemplate::compute_analytical_gradient"
void  LPTemplate::compute_analytical_gradient(PatchData &patch,
                                              Vector3D *const &grad,
                                              MsqError &err, int array_size)
{
   MsqMeshEntity* elems=patch.get_element_array(err);
   MsqVertex* vertices=patch.get_vertex_array(err);
    //Check to make sure that num_free_vert == array_size
  int num_free_vert=patch.num_free_vertices(err);
  if(array_size>=0){
    if(num_free_vert!=array_size){
      err.set_msg("Analytical Gradient passed arrays of incorrect size");
      MSQ_CHKERR(err);
    }
  }
    
  double big_f=0;
//  double total_value=0;
  double temp_value=0;
  int index=0;
//  double accum=0;
  
    //Set currentQM to be the first quality metric* in the list
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL)
    err.set_msg("LPTemplate has NULL QualityMetric pointer.");
  enum QualityMetric::EvaluationMode qm_mode=currentQM->get_evaluation_mode();
  int num_elements=patch.num_elements();
  int num_vertices=patch.num_vertices();
  int total_num=0;
  if(currentQM->get_evaluation_mode()!=QualityMetric::VERTEX)
    total_num=num_elements;
  else
    total_num=num_vertices;
  
  double *metric_values=new double[total_num];
    //If element-based, fill array with element quality metrics
  if(qm_mode!=QualityMetric::VERTEX){
    for (index=0; index<num_elements;++index){
      metric_values[index]=fabs(currentQM->evaluate_element(patch,
                                                            &elems[index],
                                                            err));
    }
  }
  else{//else vertex based metric
    for (index=0; index<num_vertices;++index){
        //evaluate metric for this vertex
      metric_values[index]=fabs(currentQM->evaluate_vertex(patch,
                                                           &vertices[index],
                                                           err));
    }
  }
  big_f=compute_function(metric_values, total_num, err);
    //TODO should this be done without pow function call
  big_f=pow(big_f,(1-pVal));
    //if the function is negated for minimization, then so is the gradient
  big_f*=get_negate_flag();
  MsqFreeVertexIndexIterator ind(&patch, err);
  ind.reset();
    //position in patch's vertex array
  int vert_count=0;
    //corresponding position in grad array
  int grad_pos=0;
    //position in elem array
  int elem_count=0;
  Vector3D grad_vec;
  while(ind.next()){
    vert_count=ind.value();
    grad[grad_pos].set(0.0,0.0,0.0);
    temp_value=0;
    if(qm_mode!=QualityMetric::VERTEX){
      temp_value=1;
        //TODO should be done only with local elements
      for (elem_count=0; elem_count<num_elements;++elem_count){
        currentQM->compute_gradient(patch, &elems[elem_count],
                                    vertices[vert_count],
                                    grad_vec,err);
        for(index=0;index<pVal-1;++index){
          temp_value*=metric_values[elem_count];
        }
        grad[grad_pos] += temp_value*grad_vec;
      }
      
    }
    else{
      err.set_msg("Vertex based metric gradients not yet implements");
    }
    grad[grad_pos]*=big_f;
    ++grad_pos;
    
  }
  delete metric_values;
}

    //
  
	
