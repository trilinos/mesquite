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
#include "MsqMessage.hpp"
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
    //making LPTemplate default to Analytical_gradient
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  set_negate_flag(qualitymetric->get_negate_flag());
}

#undef __FUNC__
#define __FUNC__ "LPTemplate::~LPTemplate"

//Michael:  need to clean up here
LPTemplate::~LPTemplate(){

}

#undef __FUNC__
#define __FUNC__ "LPTemplate::concrete_evaluate"
bool LPTemplate::concrete_evaluate(PatchData &patch, double &fval,
                                   MsqError &err){
    //Total value of objective function
  int index=0;
  bool lp_bool=true;
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
  int num_elements=patch.num_elements();
  int num_vertices=patch.num_vertices();
  int total_num=0;
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED)   
    total_num=num_elements;
  else if (currentQM->get_metric_type()==QualityMetric::VERTEX_BASED)
    total_num=num_vertices;
  else
    err.set_msg("Make sure MetricType is initialised in concrete QualityMetric constructor.");
  double *metric_values= new double[total_num];
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED)
  {
//    MsqMeshEntity* current_ent;
    for (index=0; index<num_elements;index++)
    {
      lp_bool=currentQM->evaluate_element(patch, (&elems[index]),
                                          metric_values[index], err);
        //if element is inverted, return false without calculating fval.    
      if(!lp_bool){
        delete[] metric_values;
        fval=0.0;
        return false;
      }
      metric_values[index]=fabs(metric_values[index]);
      MSQ_DEBUG_ACTION(3,{std::cout<< "      o  Quality metric value for element "
                          << index << "\t: " << metric_values[index] << "\n";});
    }
  }
  else if(currentQM->get_metric_type()==QualityMetric::VERTEX_BASED)
  {
    MsqVertex* vertices=patch.get_vertex_array(err);
    for (index=0; index<num_vertices;index++)
    {
        //evaluate metric for this vertex
      lp_bool=currentQM->evaluate_vertex(patch, (&vertices[index]),
                                         metric_values[index], err);
        //if element is inverted, return false without calculating fval.
      if(!lp_bool){
        delete[] metric_values;
        fval=0.0;
        return false;
      }
      metric_values[index]=fabs(metric_values[index]);
    }
  }
  fval=compute_function(metric_values, total_num, err);
  fval=pow(fval, 1/((double) pVal));
  delete[] metric_values;
  return true;
}

#undef __FUNC__
#define __FUNC__ "LPTemplate::compute_analytical_gradient"
/*! \fn LPTemplate::compute_analytical_gradient(PatchData &patch,Vector3D *const &grad, MsqError &err, int array_size)
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool LPTemplate::compute_analytical_gradient(PatchData &patch,
                                              Vector3D *const &grad,
                                              MsqError &err, int array_size)
{
  //Generate vertex to element connectivity if needed
  patch.generate_vertex_to_element_data();
  //vector for storing indices of vertex's connected elems
  std::vector<size_t> elem_on_vert_ind;
  std::vector<size_t> vert_on_vert_ind;
  MsqMeshEntity* elems=patch.get_element_array(err);
  MsqVertex* vertices=patch.get_vertex_array(err);
  bool lp_bool=true;
  int num_elements=patch.num_elements();
  int num_vertices=patch.num_vertices(); 
  if( num_vertices!=array_size && array_size>=0 ){
    err.set_msg("Analytical Gradient passed arrays of incorrect size.");
    MSQ_CHKERR(err); }

  double big_f=0;
  double temp_value=0;
  int index=0;
  
  //Set currentQM to be the first quality metric* in the list
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL)
    err.set_msg("LPTemplate has NULL QualityMetric pointer.");
  enum QualityMetric::MetricType qm_type=currentQM->get_metric_type();
  int total_num=0;
  if (qm_type==QualityMetric::ELEMENT_BASED)
    total_num=num_elements;
  else if (qm_type==QualityMetric::VERTEX_BASED)
    total_num=num_vertices;
  else
    err.set_msg("Make sure MetricType is initialised in concrete QualityMetric constructor.");

  // in dF = F'F
  double *metric_values=new double[total_num];
  // fill array with element quality metrics
  if(qm_type==QualityMetric::ELEMENT_BASED){
    for (index=0; index<num_elements;++index){
      lp_bool=currentQM->evaluate_element(patch, &elems[index],
                                          metric_values[index], err);
      MSQ_CHKERR(err);
        //if element is invalid
      if(!lp_bool){
        delete[] metric_values;
        return false;
      }
      
      metric_values[index]=fabs(metric_values[index]);
    }
  }
  // fill array with vertex quality metrics
  else if (qm_type==QualityMetric::VERTEX_BASED) {
    for (index=0; index<num_vertices;++index){
      //evaluate metric for this vertex
      lp_bool=currentQM->evaluate_vertex(patch, &vertices[index],
                                         metric_values[index], err);
      MSQ_CHKERR(err);
      //if patch is invalid
      if(!lp_bool){
        delete[] metric_values;
        return false;
      }
      metric_values[index] = fabs( metric_values[index] );
    }
  }
  big_f=compute_function(metric_values, total_num, err);
  if(big_f<0.0){
    err.set_msg("In LPTemplate, normed value is negative.");
  }
  big_f=pow(big_f,((1/ (double) pVal)-1));
  big_f*=get_negate_flag(); //if function negated for minimization, so is gradient.
  
  double dummy;
  //position in elem array
  size_t elem_pos=0;
  size_t vert_pos=0;
  
  Vector3D grad_vec;
  // Loops over free vertices
  for (int m=0; m<num_vertices; ++m) {
    grad[m].set(0.,0.,0.);
    temp_value=0;
    if (vertices[m].is_free_vertex()) {
      if(qm_type==QualityMetric::ELEMENT_BASED){
 
        patch.get_vertex_element_indices(m, elem_on_vert_ind,err);
        size_t ele_num_vtces = elem_on_vert_ind.size();
        MsqVertex** ele_free_vtces = new MsqVertex*[ele_num_vtces];
        elem_pos=0;
        //while(elem_pos<num_elements){
        while(!elem_on_vert_ind.empty()){
          elem_pos=(elem_on_vert_ind.back());
          elem_on_vert_ind.pop_back();
          ele_free_vtces[0] = &vertices[m];
          currentQM->compute_element_gradient(patch, &elems[elem_pos],
                                              ele_free_vtces,
                                              &grad_vec, 1, dummy, err);
          temp_value=1;
          for(index=0;index<pVal-1;++index){
            temp_value*=metric_values[elem_pos];
          }
          //if pval is odd and met val is negative
          if(metric_values[elem_pos]<0  && pVal%2 ){
            temp_value*=(-1);
          }          
          grad[m] += temp_value*grad_vec;
          //elem_pos++;
        }
        delete []ele_free_vtces;     
      }
      else{
      
        patch.get_adjacent_vertex_indices(m, vert_on_vert_ind,err);
        //For now we compute the metric for attached vertices and this
        //vertex, the above line gives us the attached vertices.  Now,
        //we must add this vertex.
        vert_on_vert_ind.push_back(m);
        size_t vert_num_vtces = vert_on_vert_ind.size();
        MsqVertex** vert_free_vtces = new MsqVertex*[vert_num_vtces];
        vert_pos=0;
        while(!vert_on_vert_ind.empty()){
          vert_pos=(vert_on_vert_ind.back());
          vert_on_vert_ind.pop_back();
          vert_free_vtces[0] = &vertices[m];
          currentQM->compute_vertex_gradient(patch, vertices[vert_pos],
                                             vert_free_vtces,
                                             &grad_vec, 1, dummy, err);
          temp_value=1;
          for(index=0;index<pVal-1;++index){
            temp_value*=metric_values[vert_pos];
          }
          //if pval is odd and met val is negative
          if(metric_values[vert_pos]<0  && pVal%2 ){
            temp_value*=(-1);
          }          
          grad[m] += temp_value*grad_vec;
        }
        delete []vert_free_vtces;     
      
      }
      grad[m]*=big_f;
      //PRINT_INFO("  gradx = %f, grady = %f, gradz = %f\n",grad[m][0],grad[m][1],grad[m][2]);
    
    }
  }
  
  delete metric_values;
  return true;
}


