/*!
  \file   GeneralizedConditionNumberQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include <vector>
#include "GeneralizedConditionNumberQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "GeneralizedConditionNumberQualityMetric::GeneralizedConditionNumberQualityMetric"

GeneralizedConditionNumberQualityMetric::GeneralizedConditionNumberQualityMetric()
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  set_name("Generalized Condition Number");
}
/*
#undef __FUNC__
#define __FUNC__ "GeneralizedConditionNumberQualityMetric::evaluate_element"

double GeneralizedConditionNumberQualityMetric::evaluate_element(PatchData *pd, int element_index,
                                                      MsqError &err)
{
  // THIS FUNCTION SHOULD EVENTUALLY BE MADE VERY EFFICIENT BY USING
  // RAW ARRAYS DIRECTLY.  RIGHT NOW IT IS A HACK TO PORT OPT-MS

  if ( pd->get_storage_mode() != PatchData::RAW_ARRAYS ) {
    err.set_msg("Need raw arrays in evaluate element function taking pd as argument\n");
    return(0.0);
  }
  if (element_index >= pd->get_element_array_size()) {
    err.set_msg("element index exceeds element array size\n");
    return(0.0);
  }

  Vector3D* coords = pd->get_coords_array(err);
  ConnectivityArrayT element_connectivity = (pd->get_connectivity_array(err))[element_index];

  switch (element_connectivity.entity_type) {
  case TRIANGLE:
    MsqNode node1(coords[element_connectivity.indices[0]][0], 
                  coords[element_connectivity.indices[0]][1], 
                  coords[element_connectivity.indices[0]][2]);
    MsqNode node2(coords[element_connectivity.indices[1]][0], 
                  coords[element_connectivity.indices[1]][1], 
                  coords[element_connectivity.indices[1]][2]);
    MsqNode node3(coords[element_connectivity.indices[2]][0], 
                  coords[element_connectivity.indices[2]][1], 
                  coords[element_connectivity.indices[2]][2]);
    MsqTri tri(&node1,&node2,&node3);
    return( this->evaluate_element(&tri, err) );
    break;
    //  default:
    //    err.set_msg("only supporting triangles in evaluate element for now\n");
    //    return(0.0);
  }
}
*/
bool GeneralizedConditionNumberQualityMetric::evaluate_element(PatchData &pd,
                                                      MsqMeshEntity *element,
                                                               double &fval,
                                                      MsqError &err)
{
  int num_sample_points;
  bool return_flag;
  std::vector<Vector3D> sample_points;
  ElementEvaluationMode eval_mode = get_element_evaluation_mode();
  element->get_sample_points(eval_mode,sample_points,err);
  std::vector<Vector3D>::iterator iter=sample_points.begin();
    // loop over sample points
  Vector3D jacobian_vectors[3];
  int num_jacobian_vectors;
  int i=0;
  num_sample_points=sample_points.size();
  double *metric_values=new double[num_sample_points];
    //Vector3D* current_sample_point;
  for(i=0;i<num_sample_points;++i){
      // compute weighted jacobian
    element->compute_weighted_jacobian(pd, (*iter),
                                       jacobian_vectors,
                                       num_jacobian_vectors, err);
      // evaluate condition number at ith sample point
      //if 2 jacobian vectors (2D elem)
    
    return_flag=compute_condition_number(pd, element, jacobian_vectors,
                                         num_jacobian_vectors,
                                         metric_values[i],err);
    if(!return_flag){
      delete metric_values;
      return false;
    }
    
    MSQ_CHKERR(err);
    ++iter;
  }// end loop over sample points
  fval=average_metrics(metric_values,num_sample_points,err);
  MSQ_CHKERR(err);
  delete metric_values;
  return true;
}

bool GeneralizedConditionNumberQualityMetric::evaluate_vertex(PatchData &/*pd*/,
                                                              MsqVertex */*vertex*/,
                                                              double &fval,
                                                              MsqError &err)
{
  err.set_msg("Condition Number's evaluate_vertex is currently being implemented");
  fval=0.0;
  return false;
  
    /*
  fval=0.0;
  size_t this_vert = pd.get_vertex_index(vert);
  std::vector<size_t> adj_elems;
  pd.get_vertex_element_indices(this_vert, adg_elems, err);
  double num_elems = adj_elems.size();
  double *metric_values=new double[ num_elems ];
  int num_jacobian_vectors;
  Vector3D sample_point;
  Vecotr3D jacobian_vectors[3];
  
  int i;
  for ( i = 0; i<num_elems; i++){
    elems[i]->get_sample_point(vertex, sample_point, err); MSQ_CHKERR(err);
    elems[i]->compute_weighted_jacobian(&current_sample_point,
                                        jacobian_vectors,
                                        num_jacobian_vectors, err);
                                        MSQ_CHKERR(err);
    metric_values[i]=compute_condition_number(jacobian_vectors,
                                              num_jacobian_vectors, err);
                                              MSQ_CHKERR(err);
  }
  total_metric=average_metrics(metric_values, num_elems, err);
  MSQ_CHKERR(err);
  delete metric_values;
    */
}
  
