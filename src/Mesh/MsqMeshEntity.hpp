/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */

/*! \file MsqMeshEntity.hpp

\author Darryl Melander
\author Thomas Leurent
\author Michael Brewer

*/

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
#include "MsqTag.hpp"

MSQ_USE(vector);
MSQ_USE(ostream);

namespace Mesquite
{
  class PatchData;

  /*!
      \class MsqMeshEntity
      \brief MsqMeshEntity is the Mesquite object that stores information about
      the elements in the mesh.

      
  */
  class MsqMeshEntity
  {
  public:
    
    MsqMeshEntity(EntityTopology type,
                  size_t vertex_indices[])
      : mType(type), mTag(0)
      {
        memcpy(vertexIndices,
               vertex_indices,
               sizeof(size_t)*vertex_count(type));
      }
    
    MsqMeshEntity()
      : mType(MIXED), mTag(0)
      {}

      //! Destructor also deletes associated tag data. 
    ~MsqMeshEntity()
      { if (mTag) { delete mTag; mTag=0; } }

     //! This operator= makes a deep copy of the tag data. 
    MsqMeshEntity& operator=(const MsqMeshEntity& rhs) { 
      mType = rhs.mType;
      memmove(vertexIndices, rhs.vertexIndices, MSQ_MAX_NUM_VERT_PER_ENT*sizeof(size_t));
      if (rhs.mTag != 0)
        mTag = new MsqTag(*(rhs.mTag));
      else
        mTag = 0;
      return *this;
    }

      //! Returns element type
    inline EntityTopology get_element_type()
      { return mType; }
    
      //! Returns the number of vertices in this element,
      //! based on its element type.
    inline size_t vertex_count() const;
    
      //! Returns the number of vertices in this element type.
    static inline size_t vertex_count(EntityTopology type);
    
      //! gets the vertices of the mesh entity
    void get_vertex_indices(vector<size_t> &vertex_list);
    void append_vertex_indices(vector<size_t> &vertex_list);
    //! Very efficient retrieval of vertices indexes 
    //! (corresponding to the PatchData vertex array).
    inline const size_t *get_vertex_index_array() const;
    inline size_t* get_modifiable_vertex_index_array();
    
      //! Sets element data
    void set_element_type(EntityTopology type)
      { mType = type; }
    void set_vertex_index(size_t vertex_in_element, size_t vertex_patch_index);
      //! Sets the vertex indices and element type in a single function call.
      //! /param indices is an array of size vertex_count(type).
    void set( EntityTopology type, const size_t *indices);
    size_t get_vertex_index(size_t vertex_in_element);
    
      //! Sets the element tag. This will overwritte an existing tag. 
    void set_tag(MsqTag* tag) {mTag = tag;}
      //! Gets the element tag.
    MsqTag* get_tag() {assert(mTag!=0); return mTag;}
    
      //fills array of Vector3D's with the jacobian vectors and the 
      //number of jacobian vectors
    void compute_weighted_jacobian(PatchData &pd, Vector3D& sample_point,
                                   Vector3D jacobian_vectors[],
                                   short &num_jacobian_vectors,
                                   MsqError &err );
    
      //!Returns a list of sample points given an evaluationmode 
    void get_sample_points(QualityMetric::ElementEvaluationMode mode,
                           vector<Vector3D> &coords,
                           MsqError &err);

    //! Returns the centroid of the element.
    void get_centroid(Vector3D& centroid, const PatchData &pd, MsqError &err) const;
    
      //!Fills a vector<size_t> with vertices connected to the given
      //!vertex through the edges of this MsqMeshEntity.
    void get_connected_vertices(size_t vertex_index,
                                vector<size_t> &vert_indices,
                                MsqError &err);
    
      //!Computes the area of the element.
      //!The returned value is always non-negative.
    double compute_unsigned_area(PatchData &pd, MsqError &err );
    
      //!Computes the volume of the element.
      //!The returned value is always non-negative.
    double compute_unsigned_volume(PatchData &pd, MsqError &err );
    
      //!Computes the signed area of the element.
    double compute_signed_area(PatchData &pd, MsqError &err );

      //!Computes the signed volume of the element.
    double compute_signed_volume(PatchData &pd, MsqError &err );
    
      //! Uses a MeshDomain call-back function to compute the normal at the corner.
    void compute_corner_normal(const size_t corner_pt, const Vector3D &corner_vec1,
                         const Vector3D &corner_vec2, Vector3D &normal,
                         PatchData &pd, MsqError &err);

      //! Compute matrices which column are the vectors issued from a corner.
      //! Stores those corner matrices in the mTag data member.  
    void compute_corner_matrices(PatchData &pd, Matrix3D A[], int num_m3d, MsqError &err );

      //! This returns a pointer to the target matrices.
      //! An error is set if they are not available. 
    TargetMatrix* get_target_matrices(size_t &num_targets, MsqError &err);

  private:
    static void get_linear_quad_jac(Vector3D *sp,
                                    Vector3D &coord0, Vector3D &coord1,
                                    Vector3D &coord2, Vector3D &coord3,
                                    Vector3D* jac);
    
    EntityTopology mType;
    size_t vertexIndices[MSQ_MAX_NUM_VERT_PER_ENT];
    MsqTag* mTag; //!< The mTag data member is a pointer, so that the memory 
                  //!< footprint stays small when no tag is used (mTag=0).
                  //!< But when a tag is pointed to by this data member, the tag
                  //!< in fact lives with the MsqMeshEntity, i.e. copies are 
                  //!< deep and the tag is deleted when the MsqMeshEntity is.
    
      // output operator for debugging.
    friend ostream& operator<<(ostream &s, const MsqMeshEntity &E);
    
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

//   inline size_t MsqMeshEntity::vertex_count(TSTT::EntityTopology type,
//                                             MsqError& err)
//   {
//     switch (type)
//     {
//       case TSTT::POINT:
//         return 1;
//       case TSTT::LINE:
//         return 2;
//       case TSTT::TRIANGLE:
//         return 3;
//       case TSTT::QUADRILATERAL:
//       case TSTT::TETRAHEDRON:
//         return 4;
//       case TSTT::PYRAMID:
//         return 5;
//       case TSTT::PRISM:
//         return 6;
//       case TSTT::SEPTAHEDRON:
//         return 7;
//       case TSTT::HEXAHEDRON:
//         return 8;
//       case TSTT::POLYGON:
//       case TSTT::POLYHEDRON:
//       default:
//         err.set_msg("Unknown element type");
//         return 0;
//     }
//   }
  
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
  
  inline void MsqMeshEntity::set(EntityTopology type, const size_t *indices)
  {
    mType = type;
    memcpy(vertexIndices, indices, sizeof(size_t) * vertex_count(type));
  }
  
  inline const size_t *MsqMeshEntity::get_vertex_index_array() const
  { return vertexIndices; }
  
  inline size_t* MsqMeshEntity::get_modifiable_vertex_index_array()
  { return vertexIndices; }
  
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

    //! \param num_targets is set to the number of corners for the element.
  inline TargetMatrix* MsqMeshEntity::get_target_matrices(size_t &num_targets, MsqError &err)
  {
    if (mTag == 0) {
      err.set_msg("no target matrix available.");
      return 0;
    }
    else {
      num_targets = vertex_count();
      return mTag->get_targets(num_targets);
    }
  }


  
  /* ***********  I/O  **************/

  inline ostream& operator<<(ostream &s, const MsqMeshEntity &E)
  {
    size_t num_vtx = E.vertex_count();
    s << "MsqMeshEntity " << &E << " with vertices ";
    for (size_t i=0; i<num_vtx; ++i)
      s << E.vertexIndices[i] << "  ";
    s << "\n";
    return s;
  }

  
} //namespace


#endif // MsqMeshEntity_hpp
