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
#define __FUNC__ "QualityMetric::compute_analytical_gradient"
/*! 
    If that function is not over-riden in the concrete class, the base
    class function makes it default to a numerical gradient.
    \param element (MsqMeshEntity*) Element on which the QualityMetric is calculated.
    \param vertex (MsqVertex*) Vertex which is considered free for purposes of computing the gradient.
    \param gradient (Vector3D &) Vector where the gradient is stored.
*/
void QualityMetric::compute_analytical_gradient(PatchData &pd,
                                                MsqMeshEntity* element,
                                                MsqVertex &vertex,
                                                Vector3D &gradient,
                                                MsqError &err)
{
  PRINT_WARNING("QualityMetric has no analytical gradient defined. ",
                "Defaulting to numerical gradient.\n");
  set_gradient_type(NUMERICAL_GRADIENT);
  compute_numerical_gradient(pd, element, vertex, gradient, err);
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_numerical_gradient"
/*!
  Numerically calculates the gradient of the QualityMetric value on
  the given element for the given free vertex.
  
  \param element (MsqMeshEntity*) Element on which the QualityMetric is
  calculated.
  \param vertex (MsqVertex*) Vertex which is considered free for purposes of
  computing the gradient.
  \param gradient (Vector3D &) Vector where the gradient is stored.
*/
void QualityMetric::compute_numerical_gradient(PatchData &pd,
                                               MsqMeshEntity *element,
                                               MsqVertex &vertex,
                                               Vector3D &gradient,
                                               MsqError &err)
{
    /*!TODO: (MICHAEL)  Try to inline this function (currenlty conflicts
      with MsqVertex.hpp).*/    
  MSQ_DEBUG_PRINT(2,"Computing Gradient\n");
  
  double delta = 10e-6;
  double func=this->evaluate_element(pd, element, err); MSQ_CHKERR(err);
  double func1=0;
  
  /* gradient in the x, y, z direction */
  for (int j=0;j<3;++j)
  {
      // perturb the coordinates of the free vertex in the j direction by delta
    vertex[j]+=delta;
      //compute the function at the perturbed point location
    func1=this->evaluate_element(pd, element, err); MSQ_CHKERR(err);
      //compute the numerical gradient
    gradient[j]=(func1-func)/delta;
      // put the coordinates back where they belong
    vertex[j] -= delta;
  }
}
