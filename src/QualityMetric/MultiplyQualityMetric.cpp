/*!
  \file   MultiplyQualityMetric.cpp
  \brief 

  \author Michael Brewer
  \date   2002-05-09
*/

#include "MultiplyQualityMetric.hpp"
#include "QualityMetric.hpp"
#include "Vector3D.hpp"
#include "MsqMessage.hpp"

using namespace Mesquite;
#undef __FUNC__
#define __FUNC__ "MultiplyQualityMetric::MultiplyQualityMetric"

MultiplyQualityMetric::MultiplyQualityMetric(QualityMetric* qm1, QualityMetric* qm2, MsqError &err){
  if(qm1 == NULL || qm2 == NULL){
    err.set_msg("MultiplyQualityMetric constructor passed NULL pointer.");
  }
  set_qmetric2(qm2);
  set_qmetric1(qm1);
  feasible=qm1->get_feasible_constraint();
  if(qm2->get_feasible_constraint())
    feasible=qm2->get_feasible_constraint();
  int n_flag=qm1->get_negate_flag();
  if(n_flag!=qm2->get_negate_flag()){
    PRINT_WARNING("MultiplyQualityMetric is being used to compose a metric that should be minimized\n with a metric that should be maximized.");
    set_negate_flag(1);
  }
  else{
    set_negate_flag(n_flag);
  }
  set_mode(qm1->get_evaluation_mode(),err);
  set_name("Composite Multiply");
}

#undef __FUNC__
#define __FUNC__ "MultiplyQualityMetric::evaluate_element"

/*! Returns qMetric1->evaluate_element(element, err) multiplied by
  qMetric2-evaluate_element(element, err)*/
double MultiplyQualityMetric::evaluate_element(PatchData& pd,
                                               MsqMeshEntity *element,
                                               MsqError &err)
{
  double total_metric;  
  total_metric=get_qmetric1()->evaluate_element(pd, element, err);
  total_metric*=get_qmetric2()->evaluate_element(pd, element, err); 
  return total_metric;
}

#undef __FUNC__
#define __FUNC__ "MultiplyQualityMetric::evaluate_node"
/*! Returns qMetric1->evaluate_node(node, err) multiplied by
  qMetric2-evaluate_node(node, err)*/
double MultiplyQualityMetric::evaluate_vertex(PatchData& pd,
                                              MsqVertex* vert,
                                              MsqError& err)
{
  double total_metric;
  total_metric=get_qmetric1()->evaluate_vertex(pd, vert, err);
  total_metric*=get_qmetric2()->evaluate_vertex(pd, vert, err);
  return total_metric;
}



