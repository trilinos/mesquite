/*!
  \file   PowerQualityMetric.cpp
  \brief 

  \author Michael Brewer
  \date   April 1, 2003
*/

#include "PowerQualityMetric.hpp"
#include "QualityMetric.hpp"
#include "Vector3D.hpp"
#include "MsqMessage.hpp"

using namespace Mesquite;
#undef __FUNC__
#define __FUNC__ "PowerQualityMetric::PowerQualityMetric"

PowerQualityMetric::PowerQualityMetric(QualityMetric* qm1,
                                       double pow_factor,
                                       MsqError &err){
  if(qm1 == NULL){
    err.set_msg("PowerQualityMetric constructor passed NULL pointer.");
  }
  set_scalealpha(pow_factor);
  set_qmetric1(qm1);
  feasible=qm1->get_feasible_constraint();
  int n_flag=qm1->get_negate_flag();
    //If the power is negative, then negate flag needs
    //to be (-n_flag).  If the power is very, very small
    //then we just set an error.
  if(fabs(pow_factor)<MSQ_MIN)
  {
    err.set_msg("PowerQualityMetric passed a double smaller than Mesquite's minimum.");
  }
  else if(pow_factor<0)
  {
    n_flag=(-n_flag);
  }
    
  set_negate_flag(n_flag);
  
  set_metric_type(qm1->get_metric_type());
  
  set_name("Composite Power");
}

#undef __FUNC__
#define __FUNC__ "PowerQualityMetric::evaluate_element"

/*! Returns qMetric1->evaluate_element(element, err) raised to the
 scaleAlpha power.  This function uses the math function pow(...).*/
bool PowerQualityMetric::evaluate_element(PatchData& pd,
                                          MsqMeshEntity *element,
                                          double &value,
                                          MsqError &err)
{
  bool valid_flag;
  double metric1;
  valid_flag=get_qmetric1()->evaluate_element(pd, element, metric1, err);
  if(!valid_flag)
    return false;
  value = pow(metric1,get_scalealpha());
  return valid_flag;
}

#undef __FUNC__
#define __FUNC__ "PowerQualityMetric::evaluate_vertex"
/*! Returns qMetric1->evaluate_vertex(...) raised to the scaleAlpha power.
  This function uses the math function pow(...).*/
bool PowerQualityMetric::evaluate_vertex(PatchData& pd,
                                            MsqVertex* vert,
                                            double &value,
                                            MsqError& err)
{
  bool valid_flag;
  double metric1;
  valid_flag=get_qmetric1()->evaluate_vertex(pd, vert, metric1, err);
  if(!valid_flag)
    return false;
  value = pow(metric1,get_scalealpha());
  return valid_flag;
}



