/*!
  \file   ConditionNumberQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include <vector>
#include "ConditionNumberQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "ConditionNumberQualityMetric::ConditionNumberQualityMetric"

ConditionNumberQualityMetric::ConditionNumberQualityMetric()
{
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  evalMode=QualityMetric::ELEMENT_VERTICES;
  set_name("Condition Number");
}

double ConditionNumberQualityMetric::evaluate_element(PatchData &pd,
                                                      MsqMeshEntity *element,
                                                      MsqError &err)
{
  double return_val=0.0;
  std::vector<size_t> v_i;
  element->get_vertex_indices(v_i);
    //only 3 temp_vec will be sent to cond-num calculator, but the
    //additional vector3D may be needed during the calculations
  Vector3D temp_vec[5];
  MsqVertex *vertices=pd.get_vertex_array(err);
  switch(element->get_element_type()){
    case TRIANGLE:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[0]];
        //make relative to equilateral
      temp_vec[1]=((2*temp_vec[2])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      return_val=condition_number_2d(temp_vec,err);
      break;
    case QUADRILATERAL:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      return_val=condition_number_2d(temp_vec,err);
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      return_val+=condition_number_2d(temp_vec,err);
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      return_val+=condition_number_2d(temp_vec,err);
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      return_val+=condition_number_2d(temp_vec,err);
      return_val/=4.0;
      break;
    case TETRAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
        //transform to equilateral tet
      temp_vec[1]=((2*temp_vec[3])-temp_vec[0])/MSQ_SQRT_THREE;
      temp_vec[2]=((3*temp_vec[4])-temp_vec[0]-temp_vec[3])/
        (MSQ_SQRT_THREE*MSQ_SQRT_TWO);
      return_val=condition_number_3d(temp_vec,err);
      break;
    case HEXAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[4]]-vertices[v_i[0]];
      return_val=condition_number_3d(temp_vec,err);
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      temp_vec[2]=vertices[v_i[5]]-vertices[v_i[1]];
      return_val+=condition_number_3d(temp_vec,err);
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      temp_vec[2]=vertices[v_i[6]]-vertices[v_i[2]];
      return_val+=condition_number_3d(temp_vec,err);
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      temp_vec[2]=vertices[v_i[7]]-vertices[v_i[3]];
      return_val+=condition_number_3d(temp_vec,err);
      temp_vec[0]=vertices[v_i[7]]-vertices[v_i[4]];
      temp_vec[1]=vertices[v_i[5]]-vertices[v_i[4]];
      temp_vec[2]=vertices[v_i[0]]-vertices[v_i[4]];
      return_val+=condition_number_3d(temp_vec,err);
      temp_vec[0]=vertices[v_i[4]]-vertices[v_i[5]];
      temp_vec[1]=vertices[v_i[6]]-vertices[v_i[5]];
      temp_vec[2]=vertices[v_i[1]]-vertices[v_i[5]];
      return_val+=condition_number_3d(temp_vec,err);
      temp_vec[0]=vertices[v_i[5]]-vertices[v_i[6]];
      temp_vec[1]=vertices[v_i[7]]-vertices[v_i[6]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[6]];
      return_val+=condition_number_3d(temp_vec,err);
      temp_vec[0]=vertices[v_i[6]]-vertices[v_i[7]];
      temp_vec[1]=vertices[v_i[4]]-vertices[v_i[7]];
      temp_vec[2]=vertices[v_i[3]]-vertices[v_i[7]];
      return_val+=condition_number_3d(temp_vec,err);
      return_val/=8.0;
      break;
    default:
      return_val=MSQ_MAX_CAP;
  }// end switch over element type
  return return_val;
}


