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
#include "PatchData.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "UntangleBetaQualityMetric::UntangleBetaQualityMetric"
/*! \fn UntangleBetaQualityMetric::UntangleBetaQualityMetric(double bet)
  \brief For untangle beta, the constructor defaults to the SUM
  averaging method, and to the ELEMENT_VERTICES evaluation mode.
*/
UntangleBetaQualityMetric::UntangleBetaQualityMetric(double bet)
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
  avgMethod=QualityMetric::RMS;
  feasible=0;
  set_name("Untangle Beta");
  set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  mBeta=bet;
}

#undef __FUNC__
#define __FUNC__ "UntangleBetaQualityMetric::evaluate_element"
/*!Evaluate the Untangle Beta value  of the MsqMeshEntity pointed to
  by 'element'.*/
bool UntangleBetaQualityMetric::evaluate_element(PatchData &pd,
                                                MsqMeshEntity *element,
                                                 double &fval,
                                                MsqError &err){
  int num_sample_points;
  std::vector<Vector3D> sample_points;
  // get sample points
  if(element==NULL)
     err.set_msg("Function passed NULL MsqMeshEntity pointer");
  ElementEvaluationMode eval_mode = get_element_evaluation_mode();
  element->get_sample_points(eval_mode, sample_points, err);
  
    // loop over sample points
  Vector3D jacobian_vectors[3];
  int num_jacobian_vectors;
  int i=0;
  num_sample_points=sample_points.size();
  double *metric_values=new double[num_sample_points];
  double temp_var;
    //vars used for 2d
  Vector3D surface_normal, cross_vec;
  size_t vertex1;
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
      vertex1=element->get_vertex_index(0);
      pd.get_surface_normal(vertex1,surface_normal,err);
      cross_vec=jacobian_vectors[0]*jacobian_vectors[1];
        //std::cout<<"\nsurface_normal "<<surface_normal;
        //std::cout<<"\cross_vec "<<cross_vec;
      temp_var=cross_vec.length();
      if(cross_vec%surface_normal<0.0){
        temp_var*=-1;
      }
      temp_var -= mBeta;
        //std::cout<<"temp_var == "<<temp_var;
      
      metric_values[i]=fabs(temp_var)-temp_var;
    }
      //if three jacobian vectors (3D elem)
    else if(num_jacobian_vectors==3)
    {
      temp_var=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
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
  fval=average_metrics(metric_values,num_sample_points,err);
  MSQ_CHKERR(err);
  delete metric_values;
  return true;
}

