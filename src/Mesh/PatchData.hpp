#ifndef MESQUITE_PATCHDATA_HPP
#define MESQUITE_PATCHDATA_HPP
/*!
  \file   PatchData.hpp
  \brief    This file contains the PatchData class and its associated mementos.


  The PatchData class provides the mesh information and functionality to Mesquite.
  The PatchDataVerticesMemento class allows the state of a PatchData object to be saved
  in order to later restore that object to its previous state.
  
  \author Thomas Leurent
  \date   2002-01-17
*/


#include <stddef.h>
#include <stdlib.h>

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"
#include "MesquiteError.hpp"

#include "MeshInterface.hpp"

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
      //! Removes all data, but capacity is unchanged.
    void clear();
      //! Removes data, frees memory used by patch
    void reset();
    
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
    size_t num_vertices() const
      { return numVertices;}
    size_t num_elements() const
      { return numElements; }
      /*! Sets the number of vertices/elements to "new_size", allocating
          space for those entities if necessary.  If "new_size" increases
          the number of vertices/elements, the data for the additional
          elements is NOT initialized, so be sure to fill the patch with
          valid data after calling these functions.
      */
    void set_num_vertices (size_t new_size);
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
    const size_t* get_vertex_to_elem_array(MsqError &err) const;
      //! Returns the start of the vertex->element offset
      //! array (v2e_o).  For vertex i, v2e_o[i] is the
      //! index into the vertex->element array (v2e) where
      //! vertex i's data begins.  So, v2e[v2e_o[i]] gives
      //! you the number of elements vertex i is attached
      //! to, and v2e[v2e_o[i]+1] gives you the index of
      //! the first element attached to vertex i.
    const size_t* get_vertex_to_elem_offset(MsqError &err) const;
    
    MsqVertex& vertex_by_index(size_t index);
    MsqMeshEntity& element_by_index(size_t index);
    size_t get_vertex_index(MsqVertex* vertex);
    size_t get_element_index(MsqMeshEntity* element);
    
      //! Get the coordinates of vertices attached to the specified element
    void get_element_vertex_coordinates(size_t elem_index,
                                        std::vector<Vector3D> &coords,
                                        MsqError &err);
      /*! Get the indices of vertices attached to the specified element. */
    void get_element_vertex_indices(size_t elem_index,
                                    std::vector<size_t> &vertex_indices,
                                    MsqError &err);
      /*! Get the indices of the elements attached to the specified vertex. */
    void get_vertex_element_indices(size_t vertex_index,
                                    std::vector<size_t> &elem_indices,
                                    MsqError &err);
    
      /*! Get the indices of vertices that are attached to vertex (given by
        vertex_index) by an element edge.
      */
    void get_adjacent_vertex_indices(size_t vertex_index,
                                     std::vector<size_t> &vert_indices,
                                     MsqError &err);
    
    
      /*!Get the indices of entities attached to entity (given by ent_ind).
        adj_ents is filled with the indices into the entity array of elements
        adjacent to the given element via an n-dimensional entity.
        Thus, if n = 0, the entities must be connected via a vertex.
        If n = 1, the entities must be connected via an edge.
        If n = 2, the entities must be connected via a two-dimensional element.
        NOTE:  if n is 2 and the elements in the entity array are
        two-dimensional, no entities should meet this criterion.
      */
    void get_adjacent_entities_via_n_dim(int n, size_t ent_ind,
                                         std::vector<size_t> &adj_ents,
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
          the normal will be set to (0,0,0)
          Check PatchData::domain_set() is not false first.
      */
    void get_surface_normal(size_t vertex_index, Vector3D &surf_norm,
                            MsqError &err);
    
      //! moves all vertices at once according to a set of directions.
      /*! \param delta_array An array of directions, ordered like the vertices in the PatchData.
        \param nb_vtx number of vertices.  Must be equal to the number of vertices in the PatchData.
        \param scale A scalar that multiplies the vectors given in delta_array.
      */
    void move_vertices(Vector3D delta_array[], size_t nb_vtx,
                       double scale, MsqError &err);
    
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
      //!Add a soft_fixed flag to all non-fixed vertices in the patch.
    void set_all_vertices_soft_fixed(MsqError &err);
    
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
      CULL_METHOD_2 = 2,   /*!< no other culling method yet. */
      CULL_METHOD_3 = 4,
      CULL_METHOD_4 = 8
    };

//     //set the geometry information (i.e, simplifiedEngine and mGeom)
//     void set_geometry_information(GeometryEngine ge,
//                                   SimplifiedGeometryEngine* msq_geom)
//     {
//       mGeom=ge;
//       simplifiedEngine=msq_geom;
//     }
//     GeometryEngine get_geometry_engine_type()
//     {
//       return mGeom;
//     }

    PatchType type() const
      { return mType; }

    //! Sets the originating meshSet. This is normally done in MeshSet::get_next_patch().
    //! This function is only for tests purposes. 
    void set_mesh_set(MeshSet* ms)
      { meshSet = ms; }
    //! Returns the originating meshSet.
    //! This function is only for tests purposes. 
    MeshSet* get_mesh_set()
      { return meshSet; }
    
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


#undef __FUNC__
#define __FUNC__ "PatchData::clear"
  inline void PatchData::clear()
  {
    numVertices = 0;
    numElements = 0;
    v2eValid = false;
    meshSet = NULL;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::reset"
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
  
#undef __FUNC__
#define __FUNC__ "PatchData::reserve_vertex_capacity"
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
        
        memcpy(new_array,
               vertexArray,
               sizeof(MsqVertex)*numVertices);
        memcpy(new_handles_array,
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
  
#undef __FUNC__
#define __FUNC__ "PatchData::reserve_element_capacity"
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
        memcpy(new_array, elementArray,
               sizeof(MsqMeshEntity)*numElements);
        memcpy(new_handles_array, elementHandlesArray,
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

  
// #undef __FUNC__
// #define __FUNC__ "PatchData::add_vertex"
//   /*! \fn PatchData::add_vertex(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh, double* vertex_coord, bool check_redundancy, MsqError &err)

//   \brief adds a vertex to the PatchData object.

//   add_vertex() returns the index of the vertex position in the vertices array.
//   If the vertex was already in the array, the return value will be the index
//   at which the vertex was already present. The function sets an error flag if
//   not enough memory is allocated in PatchData to add a vertex.

//   \param double* vertex_coords: an array of two or three double. This must be consistent
//   with the spaceDim data member.
//   \param bool check_redundancy: set to true to check if the vertex is redundant
//   in the array. If it is redundant, the vertex will not be added and
//   the function return value will be the index at which the vertex was
//   already inserted. 
//   */
//   inline int PatchData::add_vertex(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
//                                    const double vertex_coord[3],
//                                    bool check_redundancy,
//                                    MsqError &err, 
//                                    MsqVertex::FlagMaskID flag=MsqVertex::MSQ_NO_VTX_FLAG)
//   {
//     // checks if the vertex is already in array
//     if ( check_redundancy )
//       {
//         //       MsqVertex vertex(vertex_coord[0], vertex_coord[1], vertex_coord[2]);
//         Vector3D vertex(vertex_coord[0], vertex_coord[1], vertex_coord[2]);
//         for (int i = 0; i < numVertices; i++)
//           {
//             if (vertex == vertexArray[i]) 
//               return i; // vertex was already in array
//           }
//       }
    
//     // Checks that enough memory is allocated to add a vertex.
//     if (numVertices==vertexArraySize)
//       {
//         // Shouldn't we just re-allocate?
//         err.set_msg("Vertices array is already full.");
//         return -1;
//       }
    
//     // if we get here, we add the vertex to the array
//     vertexArray[numVertices].set(vertex_coord[0], vertex_coord[1], vertex_coord[2]);
//     vertexArray[numVertices].set_vertex_flag(flag);
//     vertexHandlesArray[numVertices].mesh = mh;
//     vertexHandlesArray[numVertices].entity = eh;
//     ++numVertices;
    
//     return (numVertices-1); // index of added vertex in 0-based arrays 
//   }
  
  
// #undef __FUNC__
// #define __FUNC__ "PatchData::add_vertex"
//   inline int PatchData::add_vertex(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
//                                    double x, double y, double z,
//                                    bool check_redundancy, MsqError &err, 
//                                    MsqVertex::FlagMaskID flag=MsqVertex::MSQ_NO_VTX_FLAG)
//   {
//     double coords[3];
//     coords[0] = x; coords[1] = y; coords[2] = z;
//     int index = add_vertex(mh, eh, coords, check_redundancy, err, flag); MSQ_CHKERR(err);
//     return index;
//   }
  

#undef __FUNC__
#define __FUNC__ "PatchData::get_coords_array"
  /*! \fn PatchData::get_vertex_array(MsqError &err) const 

  \brief Returns an array of all vertices in the PatchData.
  */
  inline MsqVertex* PatchData::get_vertex_array(MsqError &err) const 
  {
    if (vertexArray==0) 
      err.set_msg("\nWARNING: no vertex array defined.\n");
    return vertexArray;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::get_element_array" 
  /*! \fn PatchData::get_element_array(MsqError &err) const 

  \brief Returns the PatchData element array.
  */
  inline MsqMeshEntity* PatchData::get_element_array(MsqError &err) const
  {
    if (elementArray==0) 
      err.set_msg("\nWARNING: no element array defined.\n");
    return elementArray;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::set_vertex_coordinates"
  /*! \fn PatchData::set_vertex_coordinates(const Vector3D &coords, size_t index, MsqError &err)

  \brief set the coordinates of a vertex in the raw array
  */
  inline void PatchData::set_vertex_coordinates(const Vector3D &coords,
                                                size_t index,
                                                MsqError &err) 
  {
    if (vertexArray==0)
      err.set_msg("\nWARNING: no coordinates array defined.\n");
    
    if (index >= numVertices)
      err.set_msg("\nWARNING: index bigger than numVertices.\n");
    
    vertexArray[index] = coords;
  }
  
  
#undef __FUNC__
#define __FUNC__ "PatchData::get_vertex_to_elem_offset" 
  /*! \fn PatchData::get_vertex_to_elem_offset(MsqError &err) const 
   */  inline const size_t* PatchData::get_vertex_to_elem_offset(MsqError &err) const
  {
    if (!v2eValid || v2eOffset==NULL) 
      err.set_msg("\nWARNING: no vertex to element data available.\n");
    return v2eOffset;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::get_vertex_to_elem_array" 
  /*! \fn PatchData::get_vertex_to_elem_array(MsqError &err) const 
   */
  inline const size_t* PatchData::get_vertex_to_elem_array(MsqError &err) const
  {
    if (!v2eValid) 
      err.set_msg("\nWARNING: no vertex to element data available.\n");
    return v2E;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::vertex_by_index"
  inline MsqVertex& PatchData::vertex_by_index(size_t index)
  {
    return vertexArray[index];
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::element_by_index"
  inline MsqMeshEntity& PatchData::element_by_index(size_t index)
  {
    return elementArray[index];
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::get_vertex_index"
  /*! gets the index of a vertex in the PatchData vertex array,
    given a pointer to the vertex. */
  inline size_t PatchData::get_vertex_index(MsqVertex* vertex)
  {
    return vertex - vertexArray;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::get_element_index"
  inline size_t PatchData::get_element_index(MsqMeshEntity* element)
  {
    return element - elementArray;
  }

  
#undef __FUNC__
#define __FUNC__ "PatchData::create_vertices_memento"
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
    memcpy(memento->vertices, vertexArray, numVertices*sizeof(MsqVertex) );
    
    return memento;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::recreate_vertices_memento"
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
    memcpy(memento->vertices, vertexArray,numVertices*sizeof(MsqVertex) );
    
    memento->numVertices = numVertices;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::set_to_vertices_memento"
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
      err.set_msg("Memento may only be used to restore the PatchData "
                  "object from which it was created.");
      return;
    }
    
    if (memento->numVertices != numVertices)
    {
      err.set_msg("Unable to restore patch coordinates.  Number of "
                  "vertices in PatchData has changed.");
      return;
    }
    
      // copies the memento array into the PatchData array.
    memcpy(vertexArray, memento->vertices, numVertices*sizeof(MsqVertex) );
  }

  inline void PatchData::set_num_vertices (size_t new_size)
  {
    Mesquite::MsqError err;
    reserve_vertex_capacity(new_size, err);
    numVertices = new_size;
  }
  
  inline void PatchData::set_num_elements(size_t new_size)
  {
    Mesquite::MsqError err;
    reserve_element_capacity(new_size, err);
    numElements = new_size;
  }
  
} // namespace


#endif
