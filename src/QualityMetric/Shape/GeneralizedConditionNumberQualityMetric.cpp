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
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  evalMode=QualityMetric::ELEMENT_VERTICES;
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
double GeneralizedConditionNumberQualityMetric::evaluate_element(PatchData &pd,
                                                      MsqMeshEntity *element,
                                                      MsqError &err)
{
  int num_sample_points;
  std::vector<Vector3D> sample_points;
  element->get_sample_points(evalMode,sample_points,err);
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
    
    metric_values[i]=compute_condition_number(jacobian_vectors,
                                              num_jacobian_vectors, err);
    MSQ_CHKERR(err);
    ++iter;
  }// end loop over sample points
  double total_metric=average_metrics(metric_values,num_sample_points,err);
  MSQ_CHKERR(err);
  delete metric_values;
  return total_metric;
}

double GeneralizedConditionNumberQualityMetric::evaluate_vertex(PatchData &pd,
                                                     MsqVertex *vertex,
                                                     MsqError &err)
{
  err.set_msg("Condition Number's evaluate_vertex is currently being implemented");
  double total_metric=0.0;
    /*Commented out until we have new data structures
  MsqMeshEntity* elems = vertex->get_elements(err);
  double num_elems = vertex->get_num_adjacent_elements();
  double *metric_values=new double[ num_elems ];
  int num_jacobian_vectors;
  Vector3D sample_point;
  Vecotr3D jacobian_vectors[3];
    //!\TODO: Probably need to check the dimension here.  Depending on
    //   how we store the elements connected to the vertex.  If we are
    //   allowing 2-d elements in a 3-d mesh, then we don't want to use
    //   the 2-d elements here, unless we are doing a surface smooth,
    //   in which case we don't want the 3-d???
  
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
  return total_metric;
}
  
