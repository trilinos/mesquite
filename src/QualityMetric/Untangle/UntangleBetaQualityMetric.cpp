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
  
  double met_vals[MSQ_MAX_NUM_VERT_PER_ENT];
  fval=MSQ_MAX_CAP;
  const size_t* v_i = element->get_vertex_index_array();
  size_t e_ind = pd.get_element_index(element);
    //only 3 temp_vec will be sent to untangle calculator, but the
    //additional vector3Ds may be needed during the calculations
  Vector3D temp_vec[6];
  MsqVertex *vertices=pd.get_vertex_array(err);
  switch(element->get_element_type()){
    case TRIANGLE:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[0]];
        //make relative to equilateral
      temp_vec[1]=((2*temp_vec[2])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      untangle_function_2d(temp_vec,e_ind,pd,fval,err);
      return true;
    case QUADRILATERAL:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      untangle_function_2d(temp_vec,e_ind,pd,met_vals[0],err);
      
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      untangle_function_2d(temp_vec,e_ind,pd,met_vals[1],err);
      
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      untangle_function_2d(temp_vec,e_ind,pd,met_vals[2],err);
      
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      untangle_function_2d(temp_vec,e_ind,pd,met_vals[3],err);
      fval=average_metrics(met_vals, 4, err);
  return true;
    case TETRAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
        //transform to equilateral tet
      temp_vec[1]=((2*temp_vec[3])-temp_vec[0])/MSQ_SQRT_THREE;
      temp_vec[2]=((3*temp_vec[4])-temp_vec[0]-temp_vec[3])/
        (MSQ_SQRT_THREE*MSQ_SQRT_TWO);
      untangle_function_3d(temp_vec,fval,err);
      return true;
    case HEXAHEDRON:
        //transform to v_i[0]
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[4]]-vertices[v_i[0]];
      untangle_function_3d(temp_vec,met_vals[0],err);
      
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      temp_vec[2]=vertices[v_i[5]]-vertices[v_i[1]];
      untangle_function_3d(temp_vec,met_vals[1],err);
      
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      temp_vec[2]=vertices[v_i[6]]-vertices[v_i[2]];
      untangle_function_3d(temp_vec,met_vals[2],err);
      
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      temp_vec[2]=vertices[v_i[7]]-vertices[v_i[3]];
      untangle_function_3d(temp_vec,met_vals[3],err);
      
      temp_vec[0]=vertices[v_i[7]]-vertices[v_i[4]];
      temp_vec[1]=vertices[v_i[5]]-vertices[v_i[4]];
      temp_vec[2]=vertices[v_i[0]]-vertices[v_i[4]];
      untangle_function_3d(temp_vec,met_vals[4],err);
      
      temp_vec[0]=vertices[v_i[4]]-vertices[v_i[5]];
      temp_vec[1]=vertices[v_i[6]]-vertices[v_i[5]];
      temp_vec[2]=vertices[v_i[1]]-vertices[v_i[5]];
      untangle_function_3d(temp_vec,met_vals[5],err);
      
      temp_vec[0]=vertices[v_i[5]]-vertices[v_i[6]];
      temp_vec[1]=vertices[v_i[7]]-vertices[v_i[6]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[6]];
      untangle_function_3d(temp_vec,met_vals[6],err);
      
      temp_vec[0]=vertices[v_i[6]]-vertices[v_i[7]];
      temp_vec[1]=vertices[v_i[4]]-vertices[v_i[7]];
      temp_vec[2]=vertices[v_i[3]]-vertices[v_i[7]];
      untangle_function_3d(temp_vec,met_vals[7],err);
      fval=average_metrics(met_vals, 8, err);
      return true;
    default:
      err.set_msg("Element of incorrect type sent to UntangleBetaQualityMetric");
      return false;
  }// end switch over element type

  
}





     
