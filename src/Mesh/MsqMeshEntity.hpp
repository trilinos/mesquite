// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

//! \file MsqMeshEntity.hpp

#ifndef MSQMESHENTITY_HPP
#define MSQMESHENTITY_HPP

#include <vector>
#ifdef USE_C_PREFIX_INCLUDES
#include <cstring>
#include <cassert>
#else
#include <string.h>
#include <assert.h>
#endif

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityMetric.hpp"

namespace Mesquite
{
  /*!
      \class MsqMeshEntity
      \brief MsqMeshEntity is the Mesquite object that stores information about
      the elements in the mesh.*/
  class PatchData;
  class MsqMeshEntity
  {
  public:
    
    MsqMeshEntity(EntityTopology type,
                  size_t vertex_indices[])
        : mType(type)
      {
        memcpy(vertexIndices,
               vertex_indices,
               sizeof(size_t)*vertex_count(type));
      }
    MsqMeshEntity()
        : mType(MIXED)
      {}
    ~MsqMeshEntity()
      {}
    
      //! Returns element type
    inline EntityTopology get_element_type()
      { return mType; }
    
      //! Returns the number of vertices in this element,
      //! based on its element type.
    inline size_t vertex_count() const;
    
      //! Returns the number of vertices in this element type.
    static inline size_t vertex_count(EntityTopology type);
    static inline size_t vertex_count(TSTT::EntityTopology type,
                                      MsqError& err);
    
      //! gets the vertices of the mesh entity
    void get_vertex_indices(std::vector<size_t> &vertex_list);
    void append_vertex_indices(std::vector<size_t> &vertex_list);
      //! Sets element data
    void set_element_type(EntityTopology type)
      { mType = type; }
    void set_vertex_index(size_t vertex_in_element, size_t vertex_patch_index);
    size_t get_vertex_index(size_t vertex_in_element);
    
      //fills array of Vector3D's with the jacobian vectors and the 
      //number of jacobian vectors
    void compute_weighted_jacobian(PatchData &pd, Vector3D& sample_point,
                                   Vector3D jacobian_vectors[],
                                   int &num_jacobian_vectors,
                                   MsqError &err );
    
      //Returns a list of sample points given an evaluationmode 
    void get_sample_points(QualityMetric::ElementEvaluationMode mode,
                           std::vector<Vector3D> &coords,
                           MsqError &err);
    
  private:
    void get_linear_quad_jac(Vector3D *sp, Vector3D &coord0, Vector3D &coord1,
                             Vector3D &coord2, Vector3D &coord3,
                             Vector3D* jac);
    
    EntityTopology mType;
    size_t vertexIndices[30];
  };
  
  inline size_t MsqMeshEntity::vertex_count(EntityTopology type)
  {
    switch (type)
    {
      case TRIANGLE:
        return 3;
      case QUADRILATERAL:
      case TETRAHEDRON:
        return 4;
      case PYRAMID:
        return 5;
      case PRISM:
        return 6;
      case SEPTAHEDRON:
        return 7;
      case HEXAHEDRON:
        return 8;
      case POLYGON:
      case POLYHEDRON:
      default:
        return 0;
    }
  }

  inline size_t MsqMeshEntity::vertex_count(TSTT::EntityTopology type,
                                            MsqError& err)
  {
    switch (type)
    {
      case TSTT::POINT:
        return 1;
      case TSTT::LINE:
        return 2;
      case TSTT::TRIANGLE:
        return 3;
      case TSTT::QUADRILATERAL:
      case TSTT::TETRAHEDRON:
        return 4;
      case TSTT::PYRAMID:
        return 5;
      case TSTT::PRISM:
        return 6;
      case TSTT::SEPTAHEDRON:
        return 7;
      case TSTT::HEXAHEDRON:
        return 8;
      case TSTT::POLYGON:
      case TSTT::POLYHEDRON:
      default:
        err.set_msg("Unknown element type");
        return 0;
    }
  }
  
    // Returns the number of vertices in this type
    // of element, or 0 if a variable number.
  inline size_t MsqMeshEntity::vertex_count() const
  { return vertex_count(mType); }

  inline void MsqMeshEntity::set_vertex_index(size_t vertex_in_element,
                                              size_t vertex_patch_index)
  {
      // Make sure we're in range
    assert(vertex_in_element < vertex_count());
      // Set the index
    vertexIndices[vertex_in_element] = vertex_patch_index;
  }
  
  inline size_t MsqMeshEntity::get_vertex_index(size_t vertex_in_element)
  {
      // Make sure we're in range
    assert(vertex_in_element < vertex_count());
      // Return the index
    return vertexIndices[vertex_in_element];
  }

  inline void MsqMeshEntity::get_linear_quad_jac(Vector3D *sp,
                                                 Vector3D &coord0,
                                                 Vector3D &coord1,
                                                 Vector3D &coord2,
                                                 Vector3D &coord3,
                                                 Vector3D* jac)
  {
    jac[0]=coord1-coord0+(*sp)[1]*(coord2+coord0-coord3-coord1);
    jac[1]=coord3-coord0+(*sp)[0]*(coord2+coord0-coord3-coord1);
  }
  
    
  
} //namespace


#endif // MsqMeshEntity_hpp
