/*! \file LaplacianQualityMetric.cpp
  \author Michael Brewer
  \date 2002-05-14
  evaluates the sum of lengths of edges connected to a node
  for two- and three-diminsional elements
*/


#include "LaplacianQualityMetric.hpp"
#include <math.h>
#include <list>
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqVertex.hpp"

using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "LaplacianQualityMetric::evaluate_node"

//Michael Notes:  This code is not what we want, it just here
//mainly as an example of a node based metric

double LaplacianQualityMetric::evaluate_node(MsqVertex* node,
                                             MsqError &err)
{
    // ********** IMPORTANT **************
    //   We need to add support for the following:
    //     std::vector<size_t> adj_node_indices;
    //     pd.nodes_adjacent_to_node(node_index, adj_node_indices);
    // ***********************************
//   std::list<MsqVertex*> adj_nodes=node->adjacent_nodes(err);
  std::list<MsqVertex*> adj_nodes;
  int num_sample_points=adj_nodes.size();
  double *metric_values=new double[num_sample_points];
  Vector3D this_node=*node;
  Vector3D current_node;
  int i;
    //compute length squared of connected edges
  for(i=0;i<num_sample_points;i++)
  {
    current_node=*(adj_nodes.front());
    adj_nodes.pop_front();
    metric_values[i]=(this_node-current_node).length_squared();
  }//end loop over sample points
  
    //average the above values (by default this sums values)
  double total_metric=average_metrics(metric_values,num_sample_points,err);
  MSQ_CHKERR(err);
  delete [] metric_values;
  return total_metric;
}

