/*!
  \file   QualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-05-14
*/

#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_vertex_analytical_gradient"
/*! If that function is not over-riden in the concrete class, the base
    class function makes it default to a numerical gradient.
    \param vertex  Vertex which is considered free for purposes of computing the gradient.
    \param grad_vec Vector where the gradient is stored.
    \param metric_value Since the metric is computed, we return it. 
    \return true if the element is valid, false otherwise.
*/
bool QualityMetric::compute_vertex_analytical_gradient(PatchData &pd,
                                                       MsqVertex &vertex,
                                                       MsqVertex* vertices[],
                                                       Vector3D grad_vec[],
                                                       int num_vtx,
                                                       double &metric_value,
                                                       MsqError &err)
{
  PRINT_WARNING("QualityMetric has no analytical gradient defined. ",
                "Defaulting to numerical gradient.\n");
  set_gradient_type(NUMERICAL_GRADIENT);
  return compute_vertex_numerical_gradient(pd, vertex, vertices, grad_vec,
                                           num_vtx, metric_value, err);
}

#undef __FUNC__
#define __FUNC__ "QualityMetric::change_metric_type"
void QualityMetric::change_metric_type(MetricType /*t*/, MsqError &err)
{
  err.set_msg("This QualityMetric's MetricType can not be changed.");
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_analytical_gradient"
/*! If that function is not over-riden in the concrete class, the base
    class function makes it default to a numerical gradient.
    \param vertices base address of an array of pointers to the element vertices which
    are considered free for purposes of computing the gradient. The quality metric
    gradient relatice to each of those vertices is computed and stored in grad_vec.
    \param grad_vec base address of an array of Vector3D where the gradient is stored.
    \param num_vtx This is the size of the vertices and gradient arrays. 
    \param metric_value Since the metric is computed, we return it. 
    \return true if the element is valid, false otherwise.
*/
bool QualityMetric::compute_element_analytical_gradient(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* vertices[], Vector3D grad_vec[],
                                             int num_vtx, double &metric_value,
                                             MsqError &err)
{
  PRINT_WARNING("QualityMetric has no analytical gradient defined. ",
                "Defaulting to numerical gradient.\n");
  set_gradient_type(NUMERICAL_GRADIENT);
  return compute_element_numerical_gradient(pd, element, vertices, grad_vec, num_vtx, metric_value, err);
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_numerical_gradient"
/*!
  Numerically calculates the gradient of the QualityMetric value on
  the given element for the given free vertices.
  
  \param vertices base address of an array of pointers to the element vertices which
  are considered free for purposes of computing the gradient. The quality metric
  gradient relatice to each of those vertices is computed and stored in grad_vec.
  \param grad_vec base address of an array of Vector3D where the gradient is stored.
  \param num_vtx This is the size of the vertices and gradient arrays. 
  \param metric_value Since the metric is computed, we return it. 
  \return true if the element is valid, false otherwise.
*/
bool QualityMetric::compute_element_numerical_gradient(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* vertices[],
                                                       Vector3D grad_vec[],
                                             int num_vtx, double &metric_value,
                                             MsqError &err)
{
    /*!TODO: (MICHAEL)  Try to inline this function (currenlty conflicts
      with MsqVertex.hpp).*/    
  MSQ_DEBUG_PRINT(2,"Computing Gradient\n");
  
  bool valid=this->evaluate_element(pd, element, metric_value, err); MSQ_CHKERR(err);

  if (!valid)
    return false;
  
  double delta = 10e-6;
  double metric_value1=0;
  for (int v=0; v<num_vtx; ++v) 
  {
    /* gradient in the x, y, z direction */
    for (int j=0;j<3;++j) 
    {
      // perturb the coordinates of the free vertex in the j direction by delta
      (*vertices[v])[j]+=delta;
      //compute the function at the perturbed point location
      this->evaluate_element(pd, element,  metric_value1, err); MSQ_CHKERR(err);
      //compute the numerical gradient
      grad_vec[v][j]=(metric_value1-metric_value)/delta;
      // put the coordinates back where they belong
      (*vertices[v])[j] -= delta;
    }
  }
  return true;
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_vertex_numerical_gradient"
/*!  Numerically calculates the gradient of a vertex-based QualityMetric
  value on the given free vertex.  The metric is evaluated at MsqVertex
  'vertex', and the gradient is calculated with respect to the degrees
  of freedom associated with MsqVertices in the 'vertices' array.
*/
bool QualityMetric::compute_vertex_numerical_gradient(PatchData &pd,
                                                      MsqVertex &vertex,
                                                      MsqVertex* vertices[],
                                                      Vector3D grad_vec[],
                                                      int num_vtx,
                                                      double &metric_value,
                                                      MsqError &err)
{
   /*!TODO: (MICHAEL)  Try to inline this function (currenlty conflicts
      with MsqVertex.hpp).*/    
  MSQ_DEBUG_PRINT(2,"Computing Gradient (QualityMetric's numeric, vertex based.\n");
  
  bool valid=this->evaluate_vertex(pd, &(vertex), metric_value, err);
  MSQ_CHKERR(err);

  if (!valid)
    return false;
  
  double delta = 10e-6;
  double metric_value1=0;
  for (int v=0; v<num_vtx; ++v) 
  {
    /* gradient in the x, y, z direction */
    for (int j=0;j<3;++j) 
    {
      // perturb the coordinates of the free vertex in the j direction by delta
      (*vertices[v])[j]+=delta;
      //compute the function at the perturbed point location
      this->evaluate_vertex(pd, &(vertex),  metric_value1, err); MSQ_CHKERR(err);
      //compute the numerical gradient
      grad_vec[v][j]=(metric_value1-metric_value)/delta;
      // put the coordinates back where they belong
      (*vertices[v])[j] -= delta;
    }
  }
  return true;  
}

