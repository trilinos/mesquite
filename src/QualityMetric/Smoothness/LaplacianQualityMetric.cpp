/*! \file LaplacianQualityMetric.cpp
  \author Michael Brewer
  \date 2002-05-14
  Evaluates the sum of lengths of edges connected to a node
  for two- and three-diminsional elements
*/


#include "LaplacianQualityMetric.hpp"
#include <math.h>
#include <list>
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "LaplacianQualityMetric::evaluate_node"

//Michael Notes:  This code is not what we want, it just here
//mainly as an example of a node based metric

bool LaplacianQualityMetric::evaluate_vertex(PatchData &pd, MsqVertex* vert,
                                             double &fval, MsqError &err)
{
  fval=0.0;
  size_t this_vert = pd.get_vertex_index(vert);
  size_t other_vert;
  std::vector<size_t> adj_verts;
  Vector3D edg;
  pd.get_adjacent_vertex_indices(this_vert,adj_verts,err);
  MsqVertex* verts = pd.get_vertex_array(err);
  while(!adj_verts.empty()){
    other_vert=adj_verts.back();
    adj_verts.pop_back();
    edg[0]=verts[this_vert][0]-verts[other_vert][0];
    edg[1]=verts[this_vert][1]-verts[other_vert][1];
    edg[2]=verts[this_vert][2]-verts[other_vert][2];
    fval+=edg.length_squared();
  }   
  return true;
}

