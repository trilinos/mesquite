/*!
  \file   QualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-05-14
*/

#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
#include "PatchData.hpp"

using namespace Mesquite;
using std::cout;
using std::endl;
using std::cerr;


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
#define __FUNC__ "QualityMetric::compute_element_analytical_hessian"
/*! If that function is not over-riden in the concrete class, the base
    class function makes it default to a numerical hessian.
    \param vertices base address of an array of pointers to the element vertices which
    are considered free for purposes of computing the hessian. The quality metric
    gradient relatice to each of those vertices is computed and stored in grad_vec.
    \param grad_vec base address of an array of Vector3D where the gradient is stored.
    \param hessian base address of an array of Matrix3D where the upper part of the
    hessian is stored.
    \param num_vtx This is the size of the vertices arrays. The gradient array has
    the size of the number of vertices in the element, regardless.  
    \param metric_value Since the metric is computed, we return it. 
    \return true if the element is valid, false otherwise.
*/
bool QualityMetric::compute_element_analytical_hessian(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* vertices[], Vector3D grad_vec[],
                                             Matrix3D hessian[],
                                             int num_vtx, double &metric_value,
                                             MsqError &err)
{
  PRINT_WARNING("QualityMetric has no analytical hessian defined. ",
                "Defaulting to numerical hessian.\n");
  set_hessian_type(NUMERICAL_HESSIAN);
  return compute_element_numerical_hessian(pd, element, vertices, grad_vec,
                                           hessian, num_vtx, metric_value, err);
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_gradient_expanded"
/*!
  Note that for this function, grad_vec should be an array of size the
  number of vertices in el, not of size num_vtx.
*/
bool QualityMetric::compute_element_gradient_expanded(PatchData &pd,
                                                      MsqMeshEntity* el,
                                                      MsqVertex* vertices[],
                                                      Vector3D grad_vec[],
                                                      int num_vtx,
                                                      double &metric_value,
                                                      MsqError &err)
{
  int i;
  bool ret;
  Vector3D* grad_vec_nz = new Vector3D[num_vtx];
  ret = compute_element_gradient(pd, el, vertices, grad_vec_nz,
                                 num_vtx, metric_value, err);
  MSQ_CHKERR(err);

  std::vector<size_t> gv_i;
  gv_i.reserve(num_vtx);
  i=0;
  for (i=0; i<num_vtx; ++i) {
    gv_i.push_back( pd.get_vertex_ptr_index(vertices[i]) );
  }
     
  std::vector<size_t> ev_i;
  el->get_vertex_indices(ev_i);

  bool inc;
  std::vector<size_t>::iterator ev;
  std::vector<size_t>::iterator gv;
  for (ev=ev_i.begin(); ev!=ev_i.end(); ++ev) {
    inc = false; i=0;
    gv = gv_i.begin();
    while (gv!=gv_i.end()) {
      if (*ev == *gv) {
        inc = true;
        cout << "inc=true for ev " << *ev << "and gv " << *gv << endl; //dbg
        break;
      }
      ++gv;
    }
    if (inc == true)
      grad_vec[*ev] = grad_vec_nz[*gv];
    else
      grad_vec[*ev] = 0;
  }
  delete []grad_vec_nz;
  return ret;
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
  FUNCTION_TIMER_START(__FUNC__);
    /*!TODO: (MICHAEL)  Try to inline this function (currenlty conflicts
      with MsqVertex.hpp).*/    
  MSQ_DEBUG_PRINT(2,"Computing Numerical Gradient\n");
  
  bool valid=this->evaluate_element(pd, element, metric_value, err); MSQ_CHKERR(err);

  if (!valid)
    return false;
  double delta = 10e-6;
  int counter=0;
  double metric_value1=0;
  for (int v=0; v<num_vtx; ++v) 
  {
    /* gradient in the x, y, z direction */
    for (int j=0;j<3;++j) 
    {
        //re-initialize variables.
      valid=false;
      delta = 10e-6;
      counter=0;
        //perturb the node and calculate gradient.  The while loop is a
        //safety net to make sure the epsilon perturbation does not take
        //the element out of the feasible region.
      while(!valid && counter<10){
          // perturb the coordinates of the free vertex in the j direction
          // by delta       
        (*vertices[v])[j]+=delta;
          //compute the function at the perturbed point location
        valid=this->evaluate_element(pd, element,  metric_value1, err);
        MSQ_CHKERR(err);
          //compute the numerical gradient
        grad_vec[v][j]=(metric_value1-metric_value)/delta;
          // put the coordinates back where they belong
        (*vertices[v])[j] -= delta;
        ++counter;
        delta/=10.0;
      }
      if(counter>=10){
        err.set_msg("Perturbing vertex by delta caused an inverted element.");
      }
      
    }
  }
  FUNCTION_TIMER_END();
  return true;
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_numerical_hessian"
/*!
  Note that for this function, grad_vec should be an array of size the
  number of vertices in el, not of size num_vtx. Entries that do not correspond
  with the vertices argument array will be null.
*/
bool QualityMetric::compute_element_numerical_hessian(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* vertices[],
                                                       Vector3D grad_vec[],
                                                      Matrix3D hessian[],
                                             int num_vtx, double &metric_value,
                                             MsqError &err)
{
  MSQ_DEBUG_PRINT(2,"Computing Numerical Hessian\n");
  
  bool valid=this->compute_element_gradient_expanded(pd, element, vertices, grad_vec,
                                    num_vtx, metric_value, err); MSQ_CHKERR(err);
  
  if (!valid)
    return false;
  
  double delta = 10e-6;
  int nve = element->vertex_count();
  std::cout << "grad_vec: \n"; for (int i=0; i<nve; ++i) std::cout << grad_vec[i] << std::endl;  //dbg
  Vector3D* grad_vec1 = new Vector3D[nve];
  int v,w,i,j,s, sum_w, mat_index;
  for (v=0; v<nve; ++v) 
  {
    for (j=0;j<3;++j) 
    {
      // perturb the coordinates of the vertex v in the j direction by delta
      (*vertices[v])[j]+=delta;
      //compute the gradient at the perturbed point location
      valid = this->compute_element_gradient_expanded(pd, element, vertices, grad_vec1,
                                     num_vtx, metric_value, err); MSQ_CHKERR(err);
      assert(valid);
      //compute the numerical Hessian
      for (w=0; w<nve; ++w) {
        if (v>=w) {
          //finite difference to get some entries of the Hessian
          Vector3D fd = (grad_vec1[w]-grad_vec[w])/delta;
          // For the block at position w,v in a matrix, we need the corresponding index
          // (mat_index) in a 1D array containing only upper triangular blocks.
          sum_w = 0;
          for (s=1; s<=w; ++s) sum_w+=s;
          mat_index = w*nve+v-sum_w;
          
          for (i=0; i<3; ++i)
            hessian[mat_index][i][j] = fd[i];   
     
        }
      }
      // put the coordinates back where they belong
      (*vertices[v])[j] -= delta;
    }
  }

  delete[] grad_vec1;

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

