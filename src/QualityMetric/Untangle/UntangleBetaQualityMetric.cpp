/*! \file   UntangleBetaQualityMetric.cpp
  UntangleBeta is an untangle quality metric which can be used to evaluate
  the quality of two- or three-dimensional elements.

  \author Michael Brewer
  \date   2002-10-10
*/

#include "UntangleBetaQualityMetric.hpp"
#include <math.h>
#include <vector>
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqMeshEntity.hpp"

#include "TSTT_C.h"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "UntangleBetaQualityMetric::UntangleBetaQualityMetric"
/*! \fn UntangleBetaQualityMetric::UntangleBetaQualityMetric()
  \brief For untangle beta, the constructor defaults to the SUM
  averaging method, and to the ELEMENT_VERTICES evaluation mode.
*/
UntangleBetaQualityMetric::UntangleBetaQualityMetric()
{
  avgMethod=QualityMetric::SUM;
  evalMode=QualityMetric::ELEMENT_VERTICES;
  feasible=0;
  set_name("Untangle Beta");
  set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  mBeta=.05;
}

#undef __FUNC__
#define __FUNC__ "UntangleBetaQualityMetric::evaluate_element"
/*!Evaluate the Untangle Beta value  of the MsqMeshEntity pointed to
  by 'element'.*/
double UntangleBetaQualityMetric::evaluate_element(PatchData &pd,
                                                MsqMeshEntity *element,
                                                MsqError &err){
  int num_sample_points;
  std::vector<Vector3D> sample_points;
  // get sample points
  if(element==NULL)
     err.set_msg("Function passed NULL MsqMeshEntity pointer");
  element->get_sample_points(evalMode, sample_points, err);
  
    // loop over sample points
  Vector3D jacobian_vectors[3];
  int num_jacobian_vectors;
  int i=0;
  num_sample_points=sample_points.size();
  double *metric_values=new double[num_sample_points];
  double temp_var;
  std::vector<Vector3D>::iterator iter=sample_points.begin();
  for(i=0;i<num_sample_points;++i)
  {
      // compute weighted jacobian
    element->compute_weighted_jacobian(pd, (*iter), jacobian_vectors,
                                       num_jacobian_vectors, err);
      // evaluate mean ratio at ith sample point
      //if 2 jacobian vectors (2D elem)    
    if(num_jacobian_vectors==2)
    {
      err.set_msg("Untangle Beta currenlty only defined for 3D entities.");
    }
      //if three jacobian vectors (3D elem)
    else if(num_jacobian_vectors==3)
    {
      temp_var=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
      cout<<"\ntemp_var="<<temp_var;
      temp_var-=mBeta;
      metric_values[i]=fabs(temp_var)-temp_var;
        
    }
    else
    {
      metric_values[i]=MSQ_MAX_CAP;
    }
    ++iter;
  }// end loop over sample points
  
  MSQ_CHKERR(err);
  double total_metric=average_metrics(metric_values,num_sample_points,err);
  MSQ_CHKERR(err);
  delete metric_values;
  cout<<"\nreturning total = "<<total_metric;
  return total_metric;
}

