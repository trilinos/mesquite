/*!
  \file   AspectRatioGammaQualityMetric.cpp
  \brief  

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

/*!File AspectRatioGammaQualityMetric.cpp
  \brief evaluates the Aspect Ratio Gamma metric for two- and
  three-diminsional 
  elements
*/


double AspectRatioGammaQualityMetric::evaluate_element(PatchData &pd,
                                                       MsqMeshEntity* element,
                                                       MsqError &err)
{
  EntityTopology entity = element->get_element_type();
  double total_metric=0;
  double vol=0;
  Vector3D temp_vec(0,0,0);
  
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
        //sum of edges squared
      temp_vec=vert[1]-vert[0];
      total_metric=temp_vec.length_squared();
      temp_vec=vert[2]-vert[0];
      total_metric+=temp_vec.length_squared();
      temp_vec=vert[1]-vert[2];
      total_metric+=temp_vec.length_squared();
        //average sum of edges squared
      total_metric/=3.0;
        //normalize to equil. and div by area
        // Where does the 2.309 below come from?
      total_metric/=(vol*2.30940108);
      break;
    case TETRAHEDRON:
      total_metric=MSQ_MAX_CAP;
      break;
    default:
      total_metric=MSQ_MAX_CAP;
      cout<<"\nEntity type: "<<entity<<" not valid for Aspect Ratio Gamma\n";
      break;
  };
  
  return total_metric;  
}
