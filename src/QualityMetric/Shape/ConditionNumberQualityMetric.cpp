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
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  set_name("Condition Number");
}

bool ConditionNumberQualityMetric::evaluate_element(PatchData &pd,
                                                    MsqMeshEntity *element,
                                                    double &fval,
                                                    MsqError &err)
{
  bool return_flag;
  double met_vals[MSQ_MAX_NUM_VERT_PER_ENT];
  fval=MSQ_MAX_CAP;
  const size_t* v_i = element->get_vertex_index_array();
    //only 3 temp_vec will be sent to cond-num calculator, but the
    //additional vector3Ds may be needed during the calculations
  Vector3D temp_vec[6];
  MsqVertex *vertices=pd.get_vertex_array(err);
  switch(element->get_element_type()){
    case TRIANGLE:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[0]];
        //make relative to equilateral
      temp_vec[1]=((2*temp_vec[2])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      return_flag=condition_number_2d(temp_vec,v_i[0],pd,fval,err);
      return return_flag;
    case QUADRILATERAL:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      return_flag=condition_number_2d(temp_vec,v_i[0],pd,met_vals[0],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      return_flag=condition_number_2d(temp_vec,v_i[1],pd,met_vals[1],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      return_flag=condition_number_2d(temp_vec,v_i[2],pd,met_vals[2],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      return_flag=condition_number_2d(temp_vec,v_i[3],pd,met_vals[3],err);
      fval = average_metrics(met_vals, 4, err);
      return return_flag;
    case TETRAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
        //transform to equilateral tet
      temp_vec[1]=((2*temp_vec[3])-temp_vec[0])/MSQ_SQRT_THREE;
      temp_vec[2]=((3*temp_vec[4])-temp_vec[0]-temp_vec[3])/
        (MSQ_SQRT_THREE*MSQ_SQRT_TWO);
      return_flag=condition_number_3d(temp_vec,fval,err);
      return return_flag;
      break;
        /*
    case PYRAMID:
        //We compute the pyramid's "condition number" by averaging
        //the 4 tet's condition numbers, where the tets are created
        //by removing one of the four base vertices from the pyramid.
        //transform to origina v_i[0]
      temp_vec[3]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[0]];
        //find AW_inverse
      temp_vec[0]=temp_vec[3];
      temp_vec[1]=temp_vec[4]-temp_vec[3];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,met_vals[0],err);
      if(!return_flag)
        return return_flag;
        //transform to origina v_i[1]
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[1]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[1]];
        //find AW_inverse
      temp_vec[0]=temp_vec[3]-temp_vec[4];
      temp_vec[1]=temp_vec[3];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,met_vals[1],err);
      if(!return_flag)
        return return_flag;
        //transform to origina v_i[1]     
      temp_vec[3]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[4]=vertices[v_i[0]]-vertices[v_i[2]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[2]];
        //find AW_inverse
      temp_vec[0]=-temp_vec[3];
      temp_vec[1]=temp_vec[3]-temp_vec[4];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,met_vals[2],err);
      if(!return_flag)
        return return_flag;
        //transform to origina v_i[1]     
      temp_vec[3]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[4]=vertices[v_i[1]]-vertices[v_i[3]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[3]];
        //find AW_inverse
      temp_vec[0]=temp_vec[4]-temp_vec[3];
      temp_vec[1]=-temp_vec[3];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,met_vals[3],err);
      fval=average_metrics(met_vals, 4, err);
      if(!return_flag)
        return return_flag;
      break;
        */
    case HEXAHEDRON:
        //transform to v_i[0]
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[4]]-vertices[v_i[0]];
      return_flag=condition_number_3d(temp_vec,met_vals[0],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      temp_vec[2]=vertices[v_i[5]]-vertices[v_i[1]];
      return_flag=condition_number_3d(temp_vec,met_vals[1],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      temp_vec[2]=vertices[v_i[6]]-vertices[v_i[2]];
      return_flag=condition_number_3d(temp_vec,met_vals[2],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      temp_vec[2]=vertices[v_i[7]]-vertices[v_i[3]];
      return_flag=condition_number_3d(temp_vec,met_vals[3],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[7]]-vertices[v_i[4]];
      temp_vec[1]=vertices[v_i[5]]-vertices[v_i[4]];
      temp_vec[2]=vertices[v_i[0]]-vertices[v_i[4]];
      return_flag=condition_number_3d(temp_vec,met_vals[4],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[4]]-vertices[v_i[5]];
      temp_vec[1]=vertices[v_i[6]]-vertices[v_i[5]];
      temp_vec[2]=vertices[v_i[1]]-vertices[v_i[5]];
      return_flag=condition_number_3d(temp_vec,met_vals[5],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[5]]-vertices[v_i[6]];
      temp_vec[1]=vertices[v_i[7]]-vertices[v_i[6]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[6]];
      return_flag=condition_number_3d(temp_vec,met_vals[6],err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[6]]-vertices[v_i[7]];
      temp_vec[1]=vertices[v_i[4]]-vertices[v_i[7]];
      temp_vec[2]=vertices[v_i[3]]-vertices[v_i[7]];
      return_flag=condition_number_3d(temp_vec,met_vals[7],err);
      fval=average_metrics(met_vals, 8, err);
      return return_flag;
    default:
      fval=MSQ_MAX_CAP;
  }// end switch over element type
  return false;
}


