/*!
  \file   CornerJacobianQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include <vector>
#include "CornerJacobianQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "CornerJacobianQualityMetric::CornerJacobianQualityMetric"

CornerJacobianQualityMetric::CornerJacobianQualityMetric()
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
  avgMethod=QualityMetric::NONE;
  feasible=0;
  set_name("Corner Jacobian Volume (or Area)");
}

bool CornerJacobianQualityMetric::evaluate_element(PatchData &pd,
                                                    MsqMeshEntity *element,
                                                    double &fval,
                                                   MsqError &err)
{
  switch(element->get_element_type()){
    case TRIANGLE:
    case QUADRILATERAL:
      fval=element->compute_unsigned_area(pd,err);
      break;
    case TETRAHEDRON:
    case HEXAHEDRON:
      fval=element->compute_unsigned_volume(pd,err);
      break;
    default:
      fval=MSQ_MAX_CAP;
  }// end switch over element type
  return true;
}


