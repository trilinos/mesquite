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
#ifndef MESQUITE_PATCHDATA_HPP
#define MESQUITE_PATCHDATA_HPP
/*!
  \file   PatchData.hpp
  \brief    This file contains the PatchData class and its associated mementos.


  The PatchData class provides the mesh information and functionality to Mesquite.
  The PatchDataVerticesMemento class allows the state of a PatchData object to be saved
  in order to later restore that object to its previous state.
  
  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-01-17
*/

#include "Mesquite.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqVertex.hpp"
#include "MeshInterface.hpp"


#ifndef MSQ_USE_OLD_C_HEADERS
#  include <cstddef>
#  include <cstdlib>
#else
#  include <stddef.h>
#  include <stdlib.h>
#endif

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <map.h>
#  include <vector.h>
#else
#  include <map>
#  include <vector>
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
   class ostream;
#else
#  include <iosfwd>
#endif


namespace Mesquite
{
  class PatchDataVerticesMemento;
  class MeshSet;
//   class SimplifiedGeometryEngine;
  
  /*! \class PatchData
    Contains all the mesh information necessary for
    one iteration of the optimization algorithms over a
    local mesh patch. */
  class PatchData
  {
  public:    
      // Constructor/Destructor
    PatchData();
    ~PatchData();

  private:
      //! Doesn't allow PatchData to be copied implicitly.
      //! Mementos such as PatchDataVerticesMemento should be used when necessary. 
    PatchData(const PatchData &pd);
      //! Doesn't allow a PatchData object to be assigned to another.
      //! Mementos such as PatchDataVerticesMemento should be used when necessary. 
    PatchData& operator=(const PatchData &pd);
    
  public:

    enum ComputedInfo {
      MIN_UNSIGNED_AREA, //!< minimum volume or area out of all elements in the patch
      MAX_UNSIGNED_AREA, //!< maximum volume or area out of all elements in the patch
      MIN_EDGE_LENGTH, //!< minimum edge length in the patch
      MAX_EDGE_LENGTH, //!< maximum edge length in the patch
      MINMAX_SIGNED_DET2D, //!< minimum and maximum corner area out of all elements in the patch
      MINMAX_SIGNED_DET3D, //!< minimum and maximum corner volume out of all elements in the patch
      AVERAGE_DET3D //!< average corner determinant out of all elements in the patch
    };

    //! This function clears the patch information such as maximum volume, etc ... 
    void clear_computed_info() { computedInfos.clear(); }
    
    //! Returns the maximum volume or area out of all the elements in the patch 
    //! This information is stored in the patch and should not decrease performance
    //! when used properly. See also PatchData::clear_computed_info() .
    void get_minmax_element_unsigned_area(double& min, double& max, MsqError &err);
    
    //! Returns delta based on the minimum and maximum corner determinant over all elements in the patch
    //! This information is stored in the patch and should not decrease performance
    //! when used properly. See also PatchData::clear_computed_info() .
    double get_barrier_delta(MsqError &err); 

    //! Returns average corner determinant over all corners in the patch
    //! This information is stored in the patch and should not decrease performance
    //! when used properly. See also PatchData::clear_computed_info() .
    double get_average_Lambda_3d(MsqError &err); 

      //! Removes all data, but capacity is unchanged.
    void clear();
      //! Removes data, frees memory used by patch
    void reset();
      //! Reorders the mesh data 
    void reorder();

      /*! Ensures that at least 'min_num' vertices can be stored
          in the private arrays without requiring further allocation. */
    void reserve_vertex_capacity (size_t min_num_vertices,
                                  MsqError &err);
    /*! Ensures that at least 'min_num' elements can be stored
        in the private arrays without requiring further allocation.
        The actual number of elements in the patch is not modified.
    */
    void reserve_element_capacity(size_t min_num_elements,
                                  MsqError &err);
      //! number of vertices in the patch. 
    size_t num_vertices() const
      { return numVertices;}
      //! number of elements in the Patch.
    size_t num_elements() const
      { return numElements; }
      //! number of elements corners in the Patch. 
    size_t num_corners() const;

      /*! Sets the number of vertices to "new_size", allocating
          space for those entities if necessary.  If "new_size" increases
          the number of vertices, the data for the additional
          vertices is NOT initialized, so be sure to fill the patch with
          valid data after calling these functions.
        */
    void set_num_vertices (size_t new_size);
      //! idem set_num_vertices for elements.  
    void set_num_elements(size_t new_size);
    
      //! Returns the number of elements in the current patch who are
      //! free to move.  This is a costly function, since we have to check
      //! the flags of all vertices in the patch.
    int num_free_vertices(MsqError &err);
    
      //! Returns a pointer to the start of the vertex array.
    MsqVertex* get_vertex_array(MsqError &err) const;
      //! Returns a pointer to the start of the element array.
    MsqMeshEntity* get_element_array(MsqError &err) const;
      //! Returns the start of the vertex->element array.
      //! For each vertex in the patch, this array holds
      //! the number of elements the vertex is attached to,
      //! followed by the indices of those elements.
    const size_t* get_vertex_to_elem_array(MsqError &err);
      //! Returns the start of the vertex->element offset
      //! array (v2e_o).  For vertex i, v2e_o[i] is the
      //! index into the vertex->element array (v2e) where
      //! vertex i's data begins.  So, v2e[v2e_o[i]] gives
      //! you the number of elements vertex i is attached
      //! to, and v2e[v2e_o[i]+1] gives you the index of
      //! the first element attached to vertex i.
    const size_t* get_vertex_to_elem_offset(MsqError &err);
    
    MsqVertex& vertex_by_index(size_t index);
    MsqMeshEntity& element_by_index(size_t index);
    size_t get_vertex_index(MsqVertex* vertex);
    size_t get_element_index(MsqMeshEntity* element);
    
      //! Get the coordinates of vertices attached to the specified element
    void get_element_vertex_coordinates(size_t elem_index,
                                        msq_std::vector<Vector3D> &coords,
                                        MsqError &err);
      /*! Get the indices of vertices of specified element. !inefficient!*/
    void get_element_vertex_indices(size_t elem_index,
                                    msq_std::vector<size_t> &vertex_indices,
                                    MsqError &err);
      /*! Get the indices of the elements attached to the specified vertex. */
    void get_vertex_element_indices(size_t vertex_index,
                                    msq_std::vector<size_t> &elem_indices,
                                    MsqError &err);
    
      /*! Get the indices of vertices that are attached to vertex (given by
        vertex_index) by an element edge.
      */
    void get_adjacent_vertex_indices(size_t vertex_index,
                                     msq_std::vector<size_t> &vert_indices,
                                     MsqError &err);
    
    
      /*! \brief Get the indices of entities attached to entity 
	(given by ent_ind).
        adj_ents is filled with the indices into the entity array of elements
        adjacent to the given element via an n-dimensional entity.
        
      */
    void get_adjacent_entities_via_n_dim(int n, size_t ent_ind,
                                         msq_std::vector<size_t> &adj_ents,
                                         MsqError &err);
    
      /*! Create the arrays that store which elements are attached
        to each node.  If you know how many total vertex uses there are,
        pass it in.  Otherwise the PatchData will calculate that number.
      */
    void generate_vertex_to_element_data(size_t num_vertex_uses = 0);

    void set_vertex_coordinates(const Vector3D &coords,
                                size_t index,
                                MsqError &err);
      /*! Adjust the position of the specified vertex so that it
          lies on its constraining domain.  The actual domain constraint
          is managed by the MeshSet's MeshDomain object.
      */
    void snap_vertex_to_domain(size_t vertex_index, MsqError &err);

    /*! Returns whether a domain is associated with the MeshSet from which
        the Patch originates.
        If false, you cannot ask for a surface normal. */
    bool domain_set()
    { return domainSet; }
    
      /*! Get the normal of the surface for a given vertex.
          Normal is returned in Vector3D &surf_norm.  If the normal cannot
          be determined, or if the underlying domain is not a surface,
          the normal will be set to (0,0,0).
          Check PatchData::domain_set() is not false first.
      */
    void get_domain_normal_at_vertex(size_t vertex_index, bool normalize,
                                     Vector3D &surf_norm,
                                     MsqError &err) const;
    
      /*! Get the normal to the domain at the centroid (projected to the
          domain) of a given element.
          Normal is returned in Vector3D &surf_norm.  If the normal cannot
          be determined, or if the underlying domain is not a surface,
          the normal will be set to (0,0,0).
          Check PatchData::domain_set() is not false first.
      */
    void get_domain_normal_at_element(size_t elem_index, Vector3D &surf_norm,
                                      MsqError &err) const;

      //! Alternative signature. Same functionality.
    void get_domain_normal_at_element(MsqMeshEntity* elem_ptr,
                                      Vector3D &surf_norm, MsqError &err) const 
    { get_domain_normal_at_element(size_t(elem_ptr-elementArray), surf_norm, err); }
    
      //! Moves free vertices and then snaps the free vertices to the domain.
      /*\param dk an array of directions, ordered like the vertices in
        the PatchData.
        \param nb_vtx number of vertices.
        \param step_size a scalar that multiplies the vectors given in dk.
      */
    void move_free_vertices_constrained(Vector3D dk[], size_t nb_vtx,
                                        double step_size, MsqError &err);
    
    /*! Moves free vertices from a memento position along a certain direction 
      and then snaps the free vertices to the domain.
      \param dk an array of directions, ordered like the vertices in
      the PatchData.
      \param nb_vtx number of vertices.
      \param step_size a scalar that multiplies the vectors given in dk.
    */
    void set_free_vertices_constrained(PatchDataVerticesMemento* memento, 
                                       Vector3D dk[], size_t nb_vtx,
                                       double step_size, MsqError &err);

      //!Calculates the distance each vertex has moved from its original
      //!position as defined by the PatchDataVerticesMememnto.
    double get_max_vertex_movement_squared(PatchDataVerticesMemento* memento,
                                           MsqError &err);
    
      //! Updates the underlying mesh (the Mesquite::Mesh implementation) with
      //! new node coordinates and flag values.
    void update_mesh(MsqError &err);
    
      //!Remove the soft_fixed flag from all vertices in the patch.
    void set_all_vertices_soft_free(MsqError &err);
      //!Add a soft_fixed flag to all vertices in the patch.
    void set_all_vertices_soft_fixed(MsqError &err);
      //!Add a soft_fixed flag to all free vertices in the patch.
    void set_free_vertices_soft_fixed(MsqError &err);

      //! This function allocates an array of Tags and an array of TargetMatrix s.
      //! It then sets the TargetMatrices on the tags and sets the tags
      //! on the PatchData MsqMeshEntity s.
      //! This is the most efficient function to use from a memory allocation
      //! speed point of vue. It works for hybrid meshes.
      //! \param alloc_inv if true, also allocates space for the inverse of the target matrices.
    void allocate_target_matrices(MsqError &err, bool alloc_inv=false);
    
      //! Fills a PatchData with the elements attached to a center vertex.
      //! Note that all entities in the sub-patch are copies of the entities
      //! in 'this' patch.  As such, moving a vertex in the sub-patch
      //! won't move the corresponding vertex in the source patch.  Also,
      //! calling 'update_mesh()' on the sub-patch WILL modify the TSTT
      //! mesh, but the source patch won't see the changes.
    void get_subpatch(size_t center_vertex_index,
                      PatchData &pd_to_fill,
                      MsqError &err);
    
      //! Creates a memento that holds the current
      //! state of the PatchData coordinates. 
    PatchDataVerticesMemento* create_vertices_memento(MsqError &err);
    
      //! reinstantiates a memento to holds the current
      //! state of the PatchData coordinates. Improves memory management.
    void recreate_vertices_memento(PatchDataVerticesMemento* memento, MsqError &err);
    
    //! Restore the PatchData coordinates to the state
    //! contained in the memento.
    void set_to_vertices_memento(PatchDataVerticesMemento* memento,
                                 MsqError &err);
    
    //!  Tells MeshSet how to retrieve the mesh entities that will be stored in PatchData.
    /*!  The PatchType is set by the QualityImprover etc... and mesquite propagates
      it to the MeshSet.
    */
    enum PatchType
      {
        UNDEFINED_PATCH_TYPE,     /*!< Default.*/
        VERTICES_ON_VERTEX_PATCH, /*!< fills PatchData with the vertices connected
                                    through edges to the center vertex. */
        ELEMENTS_ON_VERTEX_PATCH, /*!< fills PatchData with the vertices connected
                                    through elements to the center vertex. */
        GLOBAL_PATCH              /*!< Fills PatchData with all elements and vertices
                                    contained in all the meshes of the MeshSet. */
      };

    /*! \enum culling_method
      Those are the culling method available to the users.
      Developpers: The values used in that enum are used by a bitset,
      so they have to be 2-based (2,4,8,16,32, ...)
    */
    enum culling_method {
      NO_BOUNDARY_VTX = 1, /*!< removes vertices on the boundary. (i.e. with a TSTT tag "boundary"). */
      NO_INTERIOR_VTX = 2,   /*!< removes vertices that are not on the boundary */
      CULL_METHOD_3 = 4,/*!< no other culling method yet. */
      CULL_METHOD_4 = 8
    };


    PatchType type() const
      { return mType; }

    //! Sets the originating meshSet. This is normally done in MeshSet::get_next_patch().
    //! This function is only for tests purposes. 
    void set_mesh_set(MeshSet* ms);
    
    //! Returns the originating meshSet.
    MeshSet* get_mesh_set()
      { return meshSet; }
    
    //! Display the coordinates and connectivity information
    msq_stdio::ostream& operator<<( msq_stdio::ostream& stream ) const;
   
   private:

    friend class MeshSet;

    MeshSet* meshSet;
    bool domainSet;
    PatchType mType;
    
      // Member data for the "local" patch
    size_t numVertices;
    size_t numElements;
    MsqVertex *vertexArray;
    Mesquite::Mesh::VertexHandle *vertexHandlesArray;
    MsqMeshEntity *elementArray;
    Mesquite::Mesh::ElementHandle *elementHandlesArray;
    
      // Links from vertices to elements
    size_t *v2E; // vertex to element array (csr data)
    size_t *v2eOffset; // vertex to element offset array (csr offset)
    size_t *subpatchIndexArray; // Used to make sub-patches fast
    
      // memory management
    size_t vertexArraySize;
    size_t elemArraySize;
    size_t v2eSize;
    size_t v2eOffsetSize;
    bool v2eValid;
    size_t subpatchIndexSize;

      // Patch Computed Information (maxs, mins, etc ... )
    msq_std::map<ComputedInfo, double>  computedInfos; 
    
//       //geometry information
//     GeometryEngine mGeom;
//     SimplifiedGeometryEngine* simplifiedEngine;
  };
  
  
  /*! \class PatchDataVerticesMemento
    \brief Contains a copy of the coordinates of a PatchData.

    Use PatchDataVerticesMemento when you want to change the coordinates
    of a PatchData object but also have the option to restore them.
    This class can only be instantiated through PatchData::create_vertices_memento().
  */
  class PatchDataVerticesMemento
  {
  public:
    ~PatchDataVerticesMemento()
    { delete[] vertices; }
  private:
    // Constructor accessible only to originator (i.e. PatchData)
    friend class PatchData;
    PatchDataVerticesMemento() 
      : originator(0), vertices(0), numVertices(0), arraySize(0)
    {}
    
    PatchData* originator; //!< PatchData whose state is kept
    MsqVertex *vertices; //!< array of vertices
    size_t numVertices;
    size_t arraySize;
  };


  inline void PatchData::clear()
  {
    numVertices = 0;
    numElements = 0;
    v2eValid = false;
    meshSet = NULL;
    computedInfos.clear();
  }
  
  inline void PatchData::reset()
  {
    clear();
    delete [] vertexArray;
    vertexArray = NULL;
    delete [] vertexHandlesArray;
    vertexHandlesArray = NULL;
    delete [] elementArray;
    elementArray = NULL;
    delete [] elementHandlesArray;
    elementHandlesArray = NULL;
    delete [] v2E;
    v2E = NULL;
    delete [] v2eOffset;
    v2eOffset = NULL;
    vertexArraySize = 0;
    elemArraySize = 0;
    v2eSize = 0;
    v2eOffsetSize = 0;
    v2eValid = false;
    delete [] subpatchIndexArray;
    subpatchIndexSize = 0;
  }
  
  /*! \fn PatchData::reserve_vertex_capacity(size_t min_num_vertices, MsqError &err)
    
  Allocates memory for vertex array
  if current memory is inadequate.
  The data for at least the first min_num_vertices is preserved,
  but pointers to MsqVertex
  objects may no longer be valid after this call.
  */
  inline void PatchData::reserve_vertex_capacity(size_t min_num_vertices,
                                                 MsqError &/*err*/)
  {
    if ( min_num_vertices > vertexArraySize
         || min_num_vertices < vertexArraySize/10)
    {
        // Allocate the new array
      MsqVertex* new_array = new MsqVertex[min_num_vertices];
      Mesh::VertexHandle* new_handles_array = new Mesh::VertexHandle[min_num_vertices];
        // Copy as much over from the old array as we can
      if (numVertices)
      {
        if (numVertices > min_num_vertices)
          numVertices = min_num_vertices;
        
        msq_stdc::memcpy(new_array,
               vertexArray,
               sizeof(MsqVertex)*numVertices);
        msq_stdc::memcpy(new_handles_array,
               vertexHandlesArray,
               sizeof(Mesquite::Mesh::VertexHandle)*numVertices);
      }
      
        // Switch to new array
      delete[] vertexArray;
      delete[] vertexHandlesArray;
      vertexArray = new_array;
      vertexHandlesArray = new_handles_array;
      vertexArraySize = min_num_vertices;
    }
  }
  
  /*! \fn PatchData::reserve_element_capacity(size_t num_elements, MsqError &err)
    allocates memory for element array
    if current memory is inapropriate.
    As much of the current data is preserved as possible,
    but pointers to MsqMeshEntity
    objects may no longer be valid after this call.
  */
  inline void PatchData::reserve_element_capacity(size_t min_num_elements,
                                                  MsqError &/*err*/)
  {
    if ( min_num_elements > elemArraySize
         || min_num_elements < elemArraySize/10)
    {
        // Allocate the new array
      MsqMeshEntity* new_array = new MsqMeshEntity[min_num_elements];
      Mesh::ElementHandle* new_handles_array = new Mesh::ElementHandle[min_num_elements];
      
        // Copy as much over from the old array as we can
      if (numElements)
      {
        if (numElements > min_num_elements)
          numElements = min_num_elements;
        msq_stdc::memcpy(new_array, elementArray,
               sizeof(MsqMeshEntity)*numElements);
        msq_stdc::memcpy(new_handles_array, elementHandlesArray,
               sizeof(Mesquite::Mesh::ElementHandle)*numElements);
      }
      
        // Switch to new array
      delete[] elementArray;
      elementArray = new_array;
      delete[] elementHandlesArray;
      elementHandlesArray = new_handles_array;
      elemArraySize = min_num_elements;
    }
  }

  /*! \fn PatchData::get_vertex_array(MsqError &err) const 

  \brief Returns an array of all vertices in the PatchData.
  */
  inline MsqVertex* PatchData::get_vertex_array(MsqError &err) const 
  {
    if (vertexArray==0) 
      MSQ_SETERR(err)( "No vertex array defined", MsqError::INVALID_STATE );
    return vertexArray;
  }
  
  /*! \fn PatchData::get_element_array(MsqError &err) const 

  \brief Returns the PatchData element array.
  */
  inline MsqMeshEntity* PatchData::get_element_array(MsqError &err) const
  {
    if (vertexArray==0) 
      MSQ_SETERR(err)( "No element array defined", MsqError::INVALID_STATE );
    return elementArray;
  }
  
  /*! \fn PatchData::set_vertex_coordinates(const Vector3D &coords, size_t index, MsqError &err)

  \brief set the coordinates of a vertex in the raw array
  */
  inline void PatchData::set_vertex_coordinates(const Vector3D &coords,
                                                size_t index,
                                                MsqError &err) 
  {
    if (vertexArray==0) {
      MSQ_SETERR(err)( "No vertex array defined", MsqError::INVALID_STATE );
      return;
    }
    
    if (index >= numVertices) {
      MSQ_SETERR(err)( "Index bigger than numVertices.", MsqError::INVALID_ARG );
      return;
    }
    
    vertexArray[index] = coords;
  }
  
  
  /*! \fn PatchData::get_vertex_to_elem_offset(MsqError &err) const 
   */  inline const size_t* PatchData::get_vertex_to_elem_offset(MsqError &/*err*/)
  {
      // Make sure we've got the data
    if (!v2eValid || !v2eOffset)
    {
      generate_vertex_to_element_data();
    }
    return v2eOffset;
  }
  
  /*! \fn PatchData::get_vertex_to_elem_array(MsqError &err) const 
   */
  inline const size_t* PatchData::get_vertex_to_elem_array(MsqError &/*err*/) 
  {
      // Make sure we've got the data
    if (!v2eValid || !v2eOffset)
    {
      generate_vertex_to_element_data();
    }
    return v2E;
  }
  
  inline MsqVertex& PatchData::vertex_by_index(size_t index)
  {
    return vertexArray[index];
  }
  
  inline MsqMeshEntity& PatchData::element_by_index(size_t index)
  {
    return elementArray[index];
  }
  
  /*! gets the index of a vertex in the PatchData vertex array,
    given a pointer to the vertex. */
  inline size_t PatchData::get_vertex_index(MsqVertex* vertex)
  {
    return vertex - vertexArray;
  }
  
  inline size_t PatchData::get_element_index(MsqMeshEntity* element)
  {
    return element - elementArray;
  }

  
  /*! \fn PatchData::create_vertices_memento(MsqError &err)
    This function instantiate PatchDataVerticesMemento object and returns a pointer to it.
    The PatchDataVerticesMemento contains the current state of the PatchData coordinates.
    It can be used to restore the same PatchData object to those coordinates.

    It is the responsibility of the caller to discard the PatchDataVerticesMemento
    when not needed any more.
  */
  inline PatchDataVerticesMemento* PatchData::create_vertices_memento(MsqError& /*err*/)
  {
    PatchDataVerticesMemento* memento = new PatchDataVerticesMemento;
    memento->originator = this;
    if (numVertices)
      memento->vertices = new MsqVertex[numVertices];
    memento->numVertices = numVertices;
    memento->arraySize = numVertices;
     
      // Copy the coordinates
    msq_stdc::memcpy(memento->vertices, vertexArray, numVertices*sizeof(MsqVertex) );
    
    return memento;
  }
  
  /*! \fn PatchData::recreate_vertices_memento(MsqError &err)
    This function reuses an existing PatchDataVerticesMemento object.
    The PatchDataVerticesMemento contains the current state of the PatchData coordinates.
    It can be used to restore the same PatchData object to those coordinates.
    
    It is the responsibility of the caller to delete the PatchDataVerticesMemento
    when it is no longer needed.
  */
  inline void PatchData::recreate_vertices_memento(PatchDataVerticesMemento* memento, 
                                                   MsqError& /*err*/)
  {
    memento->originator = this;
    
    if ( numVertices > memento->arraySize
         || numVertices < memento->arraySize/10)
    {
      delete[] memento->vertices;
        // Allocate the new array
      memento->vertices = new MsqVertex[numVertices];
      memento->arraySize = numVertices;
    }
    
      // Copy the coordinates
    msq_stdc::memcpy(memento->vertices, vertexArray,numVertices*sizeof(MsqVertex) );
    
    memento->numVertices = numVertices;
  }
  
  /*! \fn PatchData::set_to_vertices_memento(PatchDataVerticesMemento* memento, MsqError &err)
    This function restores a PatchData object coordinates to a previous state hold in
    a PatchDataVerticesMemento object (see create_vertices_memento() ).

    The function checks whether the memento originates from this particular PatchData object.
    The function does not destroy the memento object: this is the caller responsibility.
  */
  inline void PatchData::set_to_vertices_memento(PatchDataVerticesMemento* memento,
                                                 MsqError &err)
  {
    if (memento->originator != this)
    {
      MSQ_SETERR(err)("Memento may only be used to restore the PatchData "
                      "object from which it was created.",
                      MsqError::INVALID_ARG);
      return;
    }
    
    if (memento->numVertices != numVertices)
    {
      MSQ_SETERR(err)("Unable to restore patch coordinates.  Number of "
                      "vertices in PatchData has changed.",
                      MsqError::INVALID_STATE);
      return;
    }
    
      // copies the memento array into the PatchData array.
    msq_stdc::memcpy(vertexArray, memento->vertices, numVertices*sizeof(MsqVertex) );
  }

    /*! For example, a mesh composed of 3 quads has 12 corners,
         however the quads are positionned.
         This function works for hybrid meshes (like all Mesquite functions should). */
  inline size_t PatchData::num_corners() const
  {
    size_t num_corners =0;
    for (size_t i=0; i<numElements; ++i)
    {
      num_corners += elementArray[i].vertex_count();
    }
    return num_corners;
  }

  
  inline void PatchData::set_num_vertices (size_t new_size)
  {
    Mesquite::MsqError err;
    reserve_vertex_capacity(new_size, err);
    numVertices = new_size;
  }
  
  /*! Sets the number of elements to "new_size", allocating
      space for those entities if necessary.  If "new_size" increases
      the number of elements, the data for the additional
      elements is NOT initialized, so be sure to fill the patch with
      valid data after calling these functions.
  */
  inline void PatchData::set_num_elements(size_t new_size)
  {
    Mesquite::MsqError err;
    reserve_element_capacity(new_size, err);
    numElements = new_size;
  }
  
} // namespace


#endif
