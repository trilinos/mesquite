/*!
  \file   AspectRatioGammaQualityMetric.cpp
   
  \brief Evaluates the Aspect Ratio Gamma metric for two- and
  three-diminsional simplicial elements.
  \author Michael Brewer
  \date   2002-05-09
*/

#include "AspectRatioGammaQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include <list>

#include "TSTT_Base.h"

#include "MsqMeshEntity.hpp"
#include "PatchData.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "AspectRatioGammaQualityMetric::evaluate_element"
//note that we can define this metric for other element types?
//!Evaluate aspect ratio gamma on ``element''
bool AspectRatioGammaQualityMetric::evaluate_element(PatchData &pd,
                                                     MsqMeshEntity* element,
                                                     double &fval,
                                                     MsqError &err)
{
  EntityTopology entity = element->get_element_type();
  double vol=0;
  Vector3D temp_vec(0,0,0);
  bool return_val = true;
  
    //get element's nodes
  std::vector<Vector3D> vert;
  size_t elem_index = pd.get_element_index(element);
  pd.get_element_vertex_coordinates(elem_index, vert, err);
  
  switch(entity)
  {
    case TRIANGLE:
        //area
      vol=(((vert[1]-vert[0])*(vert[2]-vert[0])).length())/2.0;
      vol=fabs(vol);
      if(vol<MSQ_MIN){
        fval=MSQ_MAX_CAP;
      }
      else{
          //sum of edges squared
        temp_vec=vert[1]-vert[0];
        fval=temp_vec.length_squared();
        temp_vec=vert[2]-vert[0];
        fval+=temp_vec.length_squared();
        temp_vec=vert[1]-vert[2];
        fval+=temp_vec.length_squared();
          //average sum of edges squared
        fval/=3.0;
          //normalize to equil. and div by area
          // 2.309... is 4/sqrt(3) (inverse of the area of an equil. tri
        fval/=(vol*fourDivRootThree);
      }
      
      break;
    case TETRAHEDRON:
      vol=(vert[1]-vert[0])%((vert[2]-vert[0])*(vert[3]-vert[0]))/6.0;
        //sum of edges squared
      if(fabs(vol)<MSQ_MIN){
        fval=MSQ_MAX_CAP;
      }
      else{
        temp_vec=vert[1]-vert[0];
        fval=temp_vec.length_squared();
        temp_vec=vert[2]-vert[0];
        fval+=temp_vec.length_squared();
        temp_vec=vert[3]-vert[0];
        fval+=temp_vec.length_squared();
        temp_vec=vert[2]-vert[1];
        fval+=temp_vec.length_squared();
        temp_vec=vert[3]-vert[1];
        fval+=temp_vec.length_squared();
        temp_vec=vert[3]-vert[2];
        fval+=temp_vec.length_squared();
          //average sum of edges squared
        fval/=6.0;
        fval=sqrt(fval);
        fval*=(fval);
        fval*=(fval);
          //normalize to equil. and div by area
        fval/=(vol*twelveDivRootTwo);
      }
      break;
    default:
      fval=MSQ_MAX_CAP;
      std::cout<<"\nEntity type: "<<entity<<" not valid for Aspect Ratio Gamma\n";
      return_val = false;
      break;
  };
  
  return return_val;
}
