/*! \file LocalSizeQualityMetric.cpp
  \author Michael Brewer
  \date April 9, 2003
  Evaluates the corner volume (or areas) of the element corners
  attached to a given vertiex and then averages those values
  together.
*/


#include "LocalSizeQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MsqMessage.hpp"
#include "MsqMeshEntity.hpp"
using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "LocalSizeQualityMetric::evaluate_node"
/*!For the given vertex, vert, with connected elements, e_i for i=1...K,
  the LocalSizeQualityMetric computes the corner volumes (or areas) of
  each e_i at the corner defined by vert.  The corner volume is defined
  as the volume of the tet defined by the edges of an element which contain
  the common vertex, vert.  That volume is then diveded by the average corner
  volume of all the element corners connected to this vertex.  For
  vertices attached to pyramid elements, this metric is undefined.
*/
bool LocalSizeQualityMetric::evaluate_vertex(PatchData &pd, MsqVertex* vert,
                                             double &fval, MsqError &err)
{
  fval=0.0;
    //get the element array
  MsqMeshEntity* elems = pd.get_element_array(err);
    //conver the MsqVertex pointer into an index
  size_t this_vert = pd.get_vertex_index(vert);
    //get the vertex to element array and the offset array
  const size_t* elem_offset = pd.get_vertex_to_elem_offset(err);
  const size_t* v_to_e_array = pd.get_vertex_to_elem_array(err);
    //find the offset for this vertex
  size_t this_offset = elem_offset[this_vert];
    //get the number of elements attached to this vertex (given by the
    //first entry in the vertex to element array)
  size_t num_elems = v_to_e_array[this_offset];
    //PRINT_INFO("\nIN LOCAL SIZE CPP, num_elements = %i",num_elems);
  
  if(num_elems <= 0){
    return true;
  }
  
    //create an array to store the local metric values before averaging
    //Can we remove this dynamic allocatio?
  double* met_vals = new double[num_elems];
    //vector to hold the other verts which form a corner.
  std::vector<size_t> other_vertices;
  other_vertices.reserve(4);
  double total_val=0.0;
  size_t i=0;
    //loop over the elements attached to this vertex
  for(i=0;i<num_elems;++i){
      //get the vertices which (with this_vert) form the corner of
      //the ith element.
    elems[v_to_e_array[this_offset+i+1]].get_connected_vertices(this_vert,
                                                              other_vertices,
                                                              err);
      ////PRINT_INFO("\nINSIDE LOCAL SIZE CPP other_vertices size = %i",other_vertices.size());
    
    switch(other_vertices.size()){
        //if a surface element, compute the corner area
      case 2:
        met_vals[i] = compute_corner_area(pd, this_vert, other_vertices[0],
                                          other_vertices[1], err);
        break;
          //if a volume element, compute the corner volume 
      case 3:
        met_vals[i] = compute_corner_volume(pd, this_vert, other_vertices[0],
                                            other_vertices[1],
                                            other_vertices[2], err);
        break;
      default:
          //otherwise, there is was an error.  Either the wrong number
          //of vertices were returned fom get_connected_vertices or
          //the element does not have the correct number of edges
          //connected to this vertex (possibly a pyramid element).
        met_vals[i]=0.0;
        err.set_msg("Incorrect number of vertices returned from get_connected_vertices.");
    };
      //keep track of total so that we can compute the linear average
    total_val+=met_vals[i];
    //PRINT_INFO("\nIN LOCAL SIZE CPP, total_val = %f, i = %i",total_val,i);
      //clear the vector of other_vertices for re-use.
    other_vertices.clear();
    //PRINT_INFO("\nIN LOCAL SIZE CPP, after clean size = %f",other_vertices.size());
    
  }
    //calculate the linear average... num_elems is non-zero here.
  total_val /= (double) num_elems;
  //PRINT_INFO("\nIN LOCAL SIZE CPP, average = %f",total_val);
    //if the average is non-zero
    //divide each entry by the linear average
  if(total_val!=0){
    for(i=0;i<num_elems;++i){
      met_vals[i]/=total_val;
    }
      //calculate fval by averaging the corner values
    fval = average_metrics(met_vals, num_elems, err);
    //PRINT_INFO("\nIN LOCAL SIZE CPP, inside if statement");
  }
  //PRINT_INFO("\nIN LOCAL SIZE CPP, fval = %f",fval);
    //clean up the dynamically allocated array
  delete []met_vals;
    //always return true... the vertex position is never invalid
  return true;
  
}

