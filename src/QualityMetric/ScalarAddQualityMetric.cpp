/*!
  \file   ScalarAddQualityMetric.cpp
  \brief 

  \author Michael Brewer
  \date   April 14, 2003
*/

#include "ScalarAddQualityMetric.hpp"
#include "QualityMetric.hpp"
#include "Vector3D.hpp"
#include "MsqMessage.hpp"

using namespace Mesquite;
#undef __FUNC__
#define __FUNC__ "ScalarAddQualityMetric::ScalarAddQualityMetric"

ScalarAddQualityMetric::ScalarAddQualityMetric(QualityMetric* qm1,
                                       double scalar_double,
                                       MsqError &err){
  if(qm1 == NULL){
    err.set_msg("ScalarAddQualityMetric constructor passed NULL pointer.");
  }
  set_scalealpha(scalar_double);
  set_qmetric1(qm1);
  feasible=qm1->get_feasible_constraint();
  int n_flag=qm1->get_negate_flag();
  set_negate_flag(n_flag);
  
  set_metric_type(qm1->get_metric_type());
  
  set_name("Scalar Add Quality Metric");
}

#undef __FUNC__
#define __FUNC__ "ScalarAddQualityMetric::evaluate_element"

/*! Returns qMetric1->evaluate_element(element, err) plus a scalar value.
 */
bool ScalarAddQualityMetric::evaluate_element(PatchData& pd,
                                          MsqMeshEntity *element,
                                          double &value,
                                          MsqError &err)
{
  bool valid_flag;
  double metric1;
  valid_flag=get_qmetric1()->evaluate_element(pd, element, metric1, err);
  if(!valid_flag)
    return false;
  value = metric1+get_scalealpha();
  return valid_flag;
}

#undef __FUNC__
#define __FUNC__ "ScalarAddQualityMetric::evaluate_vertex"
/*! Returns the sum of qMetric1->evaluate_vertex(...) and a scalar value.
 */bool ScalarAddQualityMetric::evaluate_vertex(PatchData& pd,
                                            MsqVertex* vert,
                                            double &value,
                                            MsqError& err)
{
  bool valid_flag;
  double metric1;
  valid_flag=get_qmetric1()->evaluate_vertex(pd, vert, metric1, err);
  if(!valid_flag)
    return false;
  value = metric1+get_scalealpha();
  return valid_flag;
}



