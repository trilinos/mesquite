/*! \file   MeanRatioQualityMetric.cpp
  Mean Ratio is shape quality metric which can be used to evaluate
  the quality of two- or three-dimensional elements.

  \author Michael Brewer
  \date   2002-05-09
*/

#include "MeanRatioQualityMetric.hpp"
#include <math.h>
#include <vector>
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqMeshEntity.hpp"

#include "TSTT_C.h"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "MeanRatioQualityMetric::MeanRatioQualityMetric"
/*! \fn MeanRatioQualityMetric::MeanRatioQualityMetric()
  \brief For mean ratio, the constructor defaults to the LINEAR
  averaging method, and to the ELEMENT_VERTICES evaluation mode.
*/
MeanRatioQualityMetric::MeanRatioQualityMetric()
{
  avgMethod=QualityMetric::LINEAR;
  evalMode=QualityMetric::ELEMENT_VERTICES;
  feasible=1;
  set_name("Mean Ratio");
  set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  set_negate_flag(-1);
}

#undef __FUNC__
#define __FUNC__ "MeanRatioQualityMetric::evaluate_element"
//!Evaluate the Mean Ratio of the MsqMeshEntity pointed to by 'element'.
double MeanRatioQualityMetric::evaluate_element(PatchData &pd,
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
        //Michael:: Do we need fabs here?
      temp_var=fabs((jacobian_vectors[0]*jacobian_vectors[1]).length());
      metric_values[i]=jacobian_vectors[0].length_squared();
      metric_values[i]+=jacobian_vectors[1].length_squared();
      if(temp_var>=MSQ_MIN && metric_values[i]>=MSQ_MIN) //if not degenerate
      {
        metric_values[i]=2*temp_var/metric_values[i];
      }
      else
      {
        metric_values[i]=0.0;
      }
    }
    
      //if three jacobian vectors (3D elem)
    else if(num_jacobian_vectors==3)
    {
      temp_var=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
      metric_values[i]=jacobian_vectors[0].length_squared();
      metric_values[i]+=jacobian_vectors[1].length_squared();
      metric_values[i]+=jacobian_vectors[2].length_squared();
      if(temp_var>MSQ_MIN) //if not degenerate or inverted???
      {
        metric_values[i]=(3*pow(temp_var,MSQ_TWO_THIRDS))/metric_values[i];
      }
      else
      {
        metric_values[i]=0.0;
      }
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
  return total_metric;
}

