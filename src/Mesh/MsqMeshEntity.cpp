// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD: 31-Jan-03 at 13:51:14 by Thomas Leurent
//
/*! \file MsqMeshEntity.cpp

\brief This files implements all the memory management issues related
to the copy of the original TSTT (or other maybe) mesh entity handles
into Mesquite.
That copy is of course encapsulated in the MsqMeshEntity class.
  
    \author Thomas Leurent
    \date 2002-05-16  
 */

#include "MsqMeshEntity.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"

using namespace Mesquite;

#ifdef __FUNC__
#undef __FUNC__
#endif
#define __FUNC__ "MsqMeshEntity::get_vertex_indices"
//! Gets the indices of the vertices of this element.
//! The indices are only valid in the PatchData from which
//! this element was retrieved.
//! The order of the vertices is the canonical order for this
//! element's type.
void Mesquite::MsqMeshEntity::get_vertex_indices(std::vector<size_t> &vertices)
{
  vertices.clear();
  vertices.reserve(vertex_count());
  vertices.insert(vertices.end(),
                   vertexIndices,
                   vertexIndices + vertex_count());
}

#undef __FUNC__
#define __FUNC__ "MsqMeshEntity::append_vertex_indices"
//! Gets the indices of the vertices of this element.
//! The indices are only valid in the PatchData from which
//! this element was retrieved.
//! The order of the vertices is the canonical order for this
//! element's type.
//! The indices are placed appended to the end of the list.
void Mesquite::MsqMeshEntity::append_vertex_indices(std::vector<size_t> &vertex_list)
{
  vertex_list.insert(vertex_list.end(),
                   vertexIndices,
                   vertexIndices + vertex_count());
}

#undef __FUNC__
#define __FUNC__ "MsqMeshEntity::compute_weighted_jacobian"
//fills array of Vector3D's with the jacobian vectors and the 
//number of jacobian vecotors
void Mesquite::MsqMeshEntity::compute_weighted_jacobian(PatchData &pd,
                                                        Vector3D &sample_point,
                                                        Vector3D jacobian_vectors[],
                                                        int &num_jacobian_vectors,
                                                        MsqError &err )
{
  std::vector<size_t> v_v;
  get_vertex_indices(v_v);
  MsqVertex *vertices=pd.get_vertex_array(err);
  Vector3D ideal_coords[8];
  Vector3D ideal_jac[3];
  Vector3D phys_jac[3];
//  double deter;
  switch (mType)
  {
      //Note:: For the linear tri case we do not use sample pt.
    case TRIANGLE:
      jacobian_vectors[0].set(vertices[v_v[1]]-vertices[v_v[0]]);
      jacobian_vectors[1].set((2.0*vertices[v_v[2]]-vertices[v_v[0]]-
                               vertices[v_v[1]])*MSQ_SQRT_THREE_INV);
      num_jacobian_vectors=2;
      break;
      
    case QUADRILATERAL:
        /*
      ideal_coords[0].set(0,0,0);
      ideal_coords[1].set(1,0,0);
      ideal_coords[2].set(1,1,0);
      ideal_coords[3].set(0,1,0);
      get_linear_quad_jac(sample_point,ideal_coords[0],ideal_coords[1],
                          ideal_coords[2],ideal_coords[3],ideal_jac);
      get_linear_quad_jac(sample_point,vertices[v_v[0]],vertices[v_v[1]],
                          vertices[v_v[2]],vertices[v_v[3]],phys_jac);
      deter=ideal_jac[1][1]*ideal_jac[0][0]-ideal_jac[1][0]*ideal_jac[0][1];
        //cout<<"\nideal 1:: "<<ideal_jac[0];
        //cout<<"\nideal 2:: "<<ideal_jac[1];
      jacobian_vectors[0]=(phys_jac[0]*ideal_jac[1][1]-
        phys_jac[1]*ideal_jac[0][1])/deter;
      jacobian_vectors[1]=(phys_jac[1]*ideal_jac[0][0]-
        phys_jac[0]*ideal_jac[1][0])/deter;
      
      
        */          
      
      
      jacobian_vectors[0]=(vertices[v_v[1]]-vertices[v_v[0]]+sample_point[1]*
                           (vertices[v_v[2]]+vertices[v_v[0]]-vertices[v_v[3]]-
                            vertices[v_v[1]]));
      jacobian_vectors[1]=(vertices[v_v[3]]-vertices[v_v[0]]+sample_point[0]*
                           (vertices[v_v[2]]+vertices[v_v[0]]-vertices[v_v[3]]-
                            vertices[v_v[1]]));
      num_jacobian_vectors=2;
      break;
      
    case TETRAHEDRON:
      jacobian_vectors[0]=vertices[v_v[1]]-vertices[v_v[0]];
      jacobian_vectors[1]=(2.0*vertices[v_v[2]]-vertices[v_v[0]]-
                           vertices[v_v[1]])*MSQ_SQRT_THREE_INV;
      jacobian_vectors[2]=(3.0*vertices[v_v[3]]-vertices[v_v[2]]-
                           vertices[v_v[1]]-vertices[v_v[0]])*
        MSQ_SQRT_TWO_INV*MSQ_SQRT_THREE_INV;
      num_jacobian_vectors=3;
      break;
      
    case HEXAHEDRON:
      
      jacobian_vectors[0]=vertices[v_v[1]]-vertices[v_v[0]]+
                           (sample_point[1]*(vertices[v_v[2]]+vertices[v_v[0]]-
                                             vertices[v_v[3]]-vertices[v_v[1]]))+
                           (sample_point[2]*(vertices[v_v[5]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[1]]))+
                           (sample_point[1]*sample_point[2]*(vertices[v_v[6]]+
                                                             vertices[v_v[4]]+
                                                             vertices[v_v[3]]+
                                                             vertices[v_v[1]]-
                                                             vertices[v_v[7]]-
                                                             vertices[v_v[5]]-
                                                             vertices[v_v[2]]-
                                                             vertices[v_v[0]])),
      jacobian_vectors[1]=vertices[v_v[3]]-vertices[v_v[0]]+
                           (sample_point[0]*(vertices[v_v[2]]+vertices[v_v[0]]-
                                             vertices[v_v[3]]-vertices[v_v[1]]))+
                           (sample_point[2]*(vertices[v_v[7]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[3]]))+
                           (sample_point[0]*sample_point[2]*(vertices[v_v[6]]+
                                                             vertices[v_v[4]]+
                                                             vertices[v_v[3]]+
                                                             vertices[v_v[1]]-
                                                             vertices[v_v[7]]-
                                                             vertices[v_v[5]]-
                                                             vertices[v_v[2]]-
                                                             vertices[v_v[0]]));
                           
      jacobian_vectors[2]=vertices[v_v[4]]-vertices[v_v[0]]+
                           (sample_point[0]*(vertices[v_v[5]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[1]]))+
                           (sample_point[1]*(vertices[v_v[7]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[3]]))+
                           (sample_point[0]*sample_point[1]*(vertices[v_v[6]]+
                                                             vertices[v_v[4]]+
                                                             vertices[v_v[3]]+
                                                             vertices[v_v[1]]-
                                                             vertices[v_v[7]]-
                                                             vertices[v_v[5]]-
                                                             vertices[v_v[2]]-
                                                             vertices[v_v[0]]));

      num_jacobian_vectors=3;
      break;
      
    default:
      err.set_msg("Compute_weighted_jacobian not yet defined for this entity.");
      MSQ_CHKERR(err);
  } 
  
}

#undef __FUNC__
#define __FUNC__ "MsqMeshEntity::get_sample_points"
//! Appends the coordinates of the sample point to 'coords'.
/*! Places Vector3Ds holding the sample point for a given EvaluationMode
  and element type combination into a given vector of Vector3D.
  \param QualityMetric::EvaluationMode mode Specifies the type of sample
  points being used.
  \param std::vector<Vecotr3D> &coords A vector of Vector3D passed by
  reference which is used to store the sample points.
*/
void Mesquite::MsqMeshEntity::get_sample_points(QualityMetric::ElementEvaluationMode mode,
                                      std::vector<Vector3D> &coords,
                                      MsqError &err){
  switch (mType)
  {
    case TRIANGLE:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(3);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO, 0.0));
          break;
          
            //The following need to be verified
        case (QualityMetric::LINEAR_GAUSS_POINTS):
          coords.reserve(1);
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          break;
          
        case (QualityMetric::QUADRATIC_GAUSS_POINTS):
          coords.reserve(3);
          coords.push_back(Vector3D(0.5, 0.0, 0.0));	
          coords.push_back(Vector3D(0.75, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          coords.push_back(Vector3D(0.25, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          break;
          
        case (QualityMetric::CUBIC_GAUSS_POINTS):
          coords.reserve(4);
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          coords.push_back(Vector3D(0.2, 2.0* MSQ_SQRT_THREE_DIV_TWO/15.0,
                                    0.0));
          coords.push_back(Vector3D(0.8, 2.0* MSQ_SQRT_THREE_DIV_TWO/15.0,
                                    0.0));
          coords.push_back(Vector3D(0.5, 11.0* MSQ_SQRT_THREE_DIV_TWO/15.0,
                                    0.0));
          break;
        default:
            //return error saying sample points for mode not implem.    
          err.set_msg("Requested Sample Point Mode not implemented");
      }
      break;
    
    case QUADRILATERAL:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(4);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(1.0, 1.0, 0.0));
          coords.push_back(Vector3D(0.0, 1.0, 0.0));
          break;
            //THESE NEED TO BE VERIFIED
        case (QualityMetric::LINEAR_GAUSS_POINTS):
          coords.push_back(Vector3D(0.5, 0.5, 0.0));
          break;
        case (QualityMetric::QUADRATIC_GAUSS_POINTS):
        case (QualityMetric::CUBIC_GAUSS_POINTS):
        default:
            //return error saying sample points for mode not implem.
          err.set_msg("Requested Sample Point Mode not implemented");
      }
      break;

    case TETRAHEDRON:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(4);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO, 0.0));
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO/3.0,
                                    MSQ_SQRT_TWO_DIV_SQRT_THREE));
          break;
        case (QualityMetric::LINEAR_GAUSS_POINTS):

        case (QualityMetric::QUADRATIC_GAUSS_POINTS):

        case (QualityMetric::CUBIC_GAUSS_POINTS):

        default:
            //return error saying sample points for mode not implem.    
          err.set_msg("Requested Sample Point Mode not implemented");
      }   
      break;

      case HEXAHEDRON:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(8);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(1.0, 1.0, 0.0));
          coords.push_back(Vector3D(0.0, 1.0, 0.0));
          coords.push_back(Vector3D(0.0, 0.0, 1.0));	
          coords.push_back(Vector3D(1.0, 0.0, 1.0));
          coords.push_back(Vector3D(1.0, 1.0, 1.0));
          coords.push_back(Vector3D(0.0, 1.0, 1.0));
          break;
        case (QualityMetric::LINEAR_GAUSS_POINTS):

        case (QualityMetric::QUADRATIC_GAUSS_POINTS):

        case (QualityMetric::CUBIC_GAUSS_POINTS):

        default:
            //return error saying sample points for mode not implem.    
          err.set_msg("Requested Sample Point Mode not implemented");
      }   
      break;
      
    default:
        //return error saying sample points for mode not implem.
      err.set_msg("Requested Sample Point Mode not implemented");
  }
}
#undef __FUNC__
#define __FUNC__ "MsqMeshEntity::compute_unsigned_area"
/*!
  \brief Computes the area of the given element.  Returned value is
  always non-negative.  If the entity passed is not a two-dimensional
  element, an error is set.*/
double MsqMeshEntity::compute_unsigned_area(PatchData &pd, MsqError &err) {
  MsqVertex* verts=pd.get_vertex_array(err);MSQ_CHKERR(err);
  double tem=0.0;
  switch (mType)
  {
   
    case TRIANGLE:
      tem =  fabs(((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
                   (verts[vertexIndices[2]]-verts[vertexIndices[0]])).length()
                  /2.0);
      return tem;
      
    case QUADRILATERAL:
      tem = ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
             (verts[vertexIndices[3]]-verts[vertexIndices[0]])).length();
      tem += ((verts[vertexIndices[3]]-verts[vertexIndices[2]])*
              (verts[vertexIndices[1]]-verts[vertexIndices[2]])).length();
      return fabs(tem)/2.0;
      
    default:
      err.set_msg("Invalid type of element passed to compute unsigned area.");
  }
  return 0;
}
                                            
#undef __FUNC__
#define __FUNC__ "MsqMeshEntity::compute_unsigned_volume"
/*!
  \brief Computes the volume of the given element.  Returned value is
  always non-negative.  If the entity passed is not a three-dimensional
  element, an error is set.*/
double MsqMeshEntity::compute_unsigned_volume(PatchData &pd, MsqError &err) {
  Vector3D sample_point(.5,.5,.5);
  Vector3D jac_vecs[3];
  int num_jacobian_vectors=-1;
  double tem=0;
  MsqVertex *verts = pd.get_vertex_array(err);MSQ_CHKERR(err);
  switch (mType)
  {
    case TETRAHEDRON:
      tem = (verts[vertexIndices[3]]-verts[vertexIndices[0]])%
        ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
         (verts[vertexIndices[1]]-verts[vertexIndices[0]]))/6.0;
      return fabs(tem);
      
    case HEXAHEDRON:
      compute_weighted_jacobian(pd,sample_point,jac_vecs,
                                num_jacobian_vectors, err );
      return fabs(jac_vecs[2]%(jac_vecs[0]*jac_vecs[1]));
      
    default:
      err.set_msg("Invalid type of element passed to compute unsigned volume.");
  }
  return 0;
}
                                              
    

  

#undef __FUNC__
#define __FUNC__ "MsqMeshEntity::get_connected_vertices"
/*!Appends the indices (in the vertex array) of the vertices to connected
  to vertex_array[vertex_index] to the end of the vector vert_indices.
  The connected vertices are right-hand ordered as defined by the
  entity.
  
*/
void Mesquite::MsqMeshEntity::get_connected_vertices(size_t vertex_index,
                                                     std::vector<size_t>
                                                     &vert_indices,
                                                     MsqError &/*err*/)
{
    //i iterates through elem's vertices
  int i=0;
    //index is set to the index in the vertexIndices corresponding
    //to vertex_index
  int index=-1;
  
  switch (mType)
  {
    case TRIANGLE:
      while(index<0 && i<3){
        if(vertexIndices[i]==vertex_index)
          index=i;
        ++i;
      }
      if(index>=0){
        vert_indices.push_back(vertexIndices[(index+1)%3]);
        vert_indices.push_back(vertexIndices[(index+2)%3]);
      }
      
      break;
      
    case QUADRILATERAL:
      while(index<0 && i<4){
        if(vertexIndices[i]==vertex_index)
          index=i;
        ++i;
      }
      if(index>=0){
        vert_indices.push_back(vertexIndices[(index+1)%4]);
        vert_indices.push_back(vertexIndices[(index+3)%4]);
      }
          
      break;
      
    case TETRAHEDRON:
      while(index<0 && i<4){
        if(vertexIndices[i]==vertex_index)
          index=i;
        ++i;
      }
      if(index>=0){
        vert_indices.push_back(vertexIndices[(index+1)%4]);
        vert_indices.push_back(vertexIndices[(index+2)%4]);
        vert_indices.push_back(vertexIndices[(index+3)%4]);
      }
      
      break;
      
    case HEXAHEDRON:
      while(index<0 && i<8){
        if(vertexIndices[i]==vertex_index)
          index=i;
        ++i;
      }
      
      if(index>=0 && index<4)
      {
        vert_indices.push_back(vertexIndices[(index+1)%4]);
        vert_indices.push_back(vertexIndices[(index+3)%4]);
        vert_indices.push_back(vertexIndices[(index)+4]);
      }
      else if(index>=4)
      {
        vert_indices.push_back(vertexIndices[(index+3)%4]+4);
        vert_indices.push_back(vertexIndices[(index+1)%4]+4);
        vert_indices.push_back(vertexIndices[(index)-4]);
      }
          
      break;
  }
  
}


