/*! \file EdgeLengthQualityMetric.cpp
  \author Michael Brewer
  \date 2002-05-14
  Evaluates the lengths of the edges attached to the given vertex.
  By default, the averaging method is set to SUM.
*/


#include "EdgeLengthQualityMetric.hpp"
#include <math.h>
#include <list>
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "EdgeLengthQualityMetric::evaluate_node"

bool EdgeLengthQualityMetric::evaluate_vertex(PatchData &pd, MsqVertex* vert,
                                             double &fval, MsqError &err)
{
  fval=0.0;
  size_t this_vert = pd.get_vertex_index(vert);
  size_t other_vert;
  std::vector<size_t> adj_verts;
  Vector3D edg;
  pd.get_adjacent_vertex_indices(this_vert,adj_verts,err);
  int num_sample_points=adj_verts.size();
  double *metric_values=new double[num_sample_points];
  MsqVertex* verts = pd.get_vertex_array(err);
  int point_counter=0;
  while(!adj_verts.empty()){
    other_vert=adj_verts.back();
    adj_verts.pop_back();
    edg[0]=verts[this_vert][0]-verts[other_vert][0];
    edg[1]=verts[this_vert][1]-verts[other_vert][1];
    edg[2]=verts[this_vert][2]-verts[other_vert][2];    
    metric_values[point_counter]=edg.length();
    ++point_counter;
  }
  fval=average_metrics(metric_values,num_sample_points,err);
  delete metric_values;
  
  return true;
  
}

