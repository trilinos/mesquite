/*!
  \file   PatchData.hpp
  \brief  

  This file contains the PatchData class and its associated mementos.
  The PatchData Class provides the mesh information and functionality to Mesquite.
  The mementos allows to save the state of a PatchData object in order to
  later restore that object to its previous state.
  
  \author Thomas Leurent
  \date   2002-01-17
*/


#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "MeshSet.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"
#include "MesquiteError.hpp"

// this has to be down here or it screws up MeshSet.hpp
#ifndef MESQUITE_PATCHDATA_HPP
#define MESQUITE_PATCHDATA_HPP

namespace Mesquite
{
  class PatchDataCoordsMemento;
  
    //! \class PatchData
    //!       Contains all the mesh information necessary for
    //!       one iteration of the optimization algorythms over the
    //!       local mesh patch.
  class PatchData
  {
  public:    
      // Constructor/Destructor
    PatchData();
    ~PatchData();
    
    //! Removes all data, but capacity is unchanged.
    void clear();
      //! Removes data, frees memory used by patch
    void reset();
    
      //! ensures that at least 'min_num' entities can be stored
      //! in the private arrays without requiring further allocation.
    void reserve_vertex_capacity (size_t min_num_vertices,
                                  MsqError &err);
    void reserve_element_capacity(size_t min_num_elements,
                                  MsqError &err);
    
    void set_num_free_vertices(int f_v)
      { numFreeVertices = f_v; }
    void set_num_free_elements(int f_e)
      { numFreeElements = f_e; }
    
    int add_vertex(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
                   double *vertex_coords,
                   bool check_redundancy,
                   MsqError &err);
    void add_element(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
                     int* vertex_indices,
                     EntityTopology topo,
                     MsqError &err);
    void add_triangle(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
                      size_t index_vertex1,
                      size_t index_vertex2,
                      size_t index_vertex3,
                      MsqError &err);
    
    int num_vertices() const
      { return numVertices; }
    int num_elements() const
      { return numElements; } 
    int num_free_vertices() const
      { return numFreeVertices; }
    int num_free_elements() const
      { return numFreeElements; }
    
      // Get the whole array at a time
    MsqVertex* get_vertex_array(MsqError &err) const;
    MsqMeshEntity* get_element_array(MsqError &err) const;
    size_t* get_vertex_to_elem_offset(MsqError &err) const;
    size_t* get_vertex_to_elem_array(MsqError &err) const;
    
    MsqVertex& vertex_by_index(size_t index);
    MsqMeshEntity& element_by_index(size_t index);
    size_t get_vertex_index(MsqVertex* vertex);
    size_t get_element_index(MsqMeshEntity* element);
    
      // Get the coordinates of vertices attached to the
      // specified element
    void get_element_vertex_coordinates(size_t elem_index,
                                        std::vector<Vector3D> &coords,
                                        MsqError &err);
    void get_element_vertex_indices(size_t elem_index,
                                    std::vector<size_t> &vertex_indices,
                                    MsqError &err);
    void get_vertex_element_indices(size_t vertex_index,
                                    std::vector<size_t> &elem_indices);
    
    void set_vertex_coordinates(const Vector3D &coords,
                                size_t index,
                                MsqError &err);
    
      //! advanced use functions 
//     double get_smallest_edge_length(MsqError &err);
     void move_free_vertices(Vector3D dk[], int nb_vtx,
                             double step_size, MsqError &err);
    
    // Updates the TSTT mesh with any changes made to the PatchData
    void update_mesh(MsqError &err);

    //! Creates a memento that holds the current
      //! state of the PatchData coordinates. 
    PatchDataCoordsMemento* create_coords_memento(MsqError &err);
    
      //! Restore the PatchData coordinates to the state
      //! contained in the memento.
    void set_to_coords_memento(PatchDataCoordsMemento* memento,
                               MsqError &err);
    
  private:

    struct EntityEntry
    {
      TSTT::Mesh_Handle mesh;
      TSTT::Entity_Handle entity;
      EntityEntry( TSTT::Mesh_Handle m, TSTT::Entity_Handle e ) :
        mesh(m), entity(e)
        {}
      EntityEntry( ) :
        mesh(0), entity(0)
      {}
    };
    
      // Don't allow PatchData to be copied implicitly
    PatchData(const PatchData &A);
      // copy function used in copy constructor and assignment
    void copy(const PatchData &A);
    
    friend class MeshSet;
    
      // Member data for the "local" patch
    int numVertices;
    int numElements;
    int numFreeVertices;
    int numFreeElements;
    MsqVertex *vertexArray;
    EntityEntry *vertexHandlesArray;
    MsqMeshEntity *elementArray;
    EntityEntry *elementHandlesArray;
    
      // memory management utilities
    int vertexArraySize;
    int elemArraySize;

      // Links from vertices to elements
    size_t *elemsInVertex;
    size_t *vertexToElemOffset;
  };


  /*! \class PatchDataCoordsMemento
      \brief Contains a copy of the coordinates of a PatchData.

      Use PatchDataCoordsMemento when you want to change the coordinates
      of a PatchData object but also have the option to restore them.
      This class can only be instantiated through PatchData::create_coords_memento().
    */
  class PatchDataCoordsMemento
  {
  public:
    ~PatchDataCoordsMemento()
      { delete[] coords; }
  private:
      // Constructor accessible only to originator (i.e. PatchData)
    friend class PatchData;
    PatchDataCoordsMemento() 
        : originator(0), coords(0), numVertices(0)
      {}
    
    PatchData* originator; // PatchData whose state is kept
    Vector3D *coords; // array of coordinates triplet
    int numVertices;
  };
  
#undef __FUNC__
#define __FUNC__ "PatchData::clear"
  inline void PatchData::clear()
  {
    numVertices = 0;
    numElements = 0;
    numFreeVertices = 0;
    numFreeElements = 0;
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
    vertexArraySize = 0;
    elemArraySize = 0;
    delete [] elemsInVertex;
    elemsInVertex = NULL;
    delete [] vertexToElemOffset;
    vertexToElemOffset = NULL;
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
                                                 MsqError &err)
  {
    if ( min_num_vertices > vertexArraySize
         || min_num_vertices < vertexArraySize/10)
    {
        // Allocate the new array
      MsqVertex* new_array = new MsqVertex[min_num_vertices];
      EntityEntry* new_handles_array = new EntityEntry[min_num_vertices];
        // Copy as much over from the old array as we can
      if (numVertices)
      {
        if (numVertices > min_num_vertices)
          numVertices = min_num_vertices;
        memcpy(new_array, vertexArray, sizeof(MsqVertex)*numVertices);
        memcpy(new_handles_array, vertexHandlesArray, sizeof(EntityEntry)*numVertices);
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
                                                  MsqError &err)
  {
    if ( min_num_elements > elemArraySize
         || min_num_elements < elemArraySize/10)
    {
        // Allocate the new array
      MsqMeshEntity* new_array = new MsqMeshEntity[min_num_elements];
      EntityEntry* new_handles_array = new EntityEntry[min_num_elements];
      
        // Copy as much over from the old array as we can
      if (numElements)
      {
        if (numElements > min_num_elements)
          numElements = min_num_elements;
        memcpy(new_array, elementArray,
               sizeof(MsqMeshEntity)*numElements);
        memcpy(new_handles_array, elementHandlesArray,
               sizeof(EntityEntry)*numElements);
      }
      
        // Switch to new array
      delete[] elementArray;
      elementArray = new_array;
      delete[] elementHandlesArray;
      elementHandlesArray = new_handles_array;
      elemArraySize = min_num_elements;
    }
  }

#undef __FUNC__
#define __FUNC__ "PatchData::add_vertex"
  /*! \fn PatchData::add_vertex(double* vertex_coord, bool check_redundancy, MsqError &err)

  \brief adds a vertex to the PatchData object.

  add_vertex() returns the index of the vertex position in the vertices array.
  If the vertex was already in the array, the return value will be the index
  at which the vertex was already present. The function sets an error flag if
  not enough memory is allocated in PatchData to add a vertex.

  \param double* vertex_coords: an array of two or three double. This must be consistent
         with the spaceDim data member.
  \param bool check_redundancy: set to true to check if the vertex is redundant
         in the array. If it is redundant, the vertex will not be added and
         the function return value will be the index at which the vertex was
         already inserted. 
  */
  inline int PatchData::add_vertex(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
                                   double* vertex_coord, bool check_redundancy,
                                   MsqError &err)
  {
    MsqVertex vertex(vertex_coord[0], vertex_coord[1], vertex_coord[2]);
    
      // checks if the vertex is already in an array
    if ( check_redundancy )
    {
      for (int i = 0; i < numVertices; i++)
      {
        if (vertex == vertexArray[i]) 
          return i; // vertex was already in array
      }
    }
    
      // Checks that enough memory is allocated to add a vertex.
    if (numVertices==vertexArraySize)
    {
        // Shouldn't we just re-allocate?
      err.set_msg("Vertices array is already full.");
      return -1;
    }
    
      // if we get here, we add the vertex to the array
    vertexArray[numVertices] = vertex;
    vertexHandlesArray[numVertices].mesh = mh;
    vertexHandlesArray[numVertices].entity = eh;
    ++numVertices;
    
    return (numVertices-1); // index of added vertex in 0-based arrays 
  }
  
  
#undef __FUNC__
#define __FUNC__ "PatchData::get_coords_array"
  /*! PatchData::get_coords_array(MsqError &err)

      \brief return the PatchData vertices information as a C array of doubles.

      This can only be used if the storage mode has been set to RAW_ARRAYS.
      See set_storage_mode().
    */
  inline MsqVertex* PatchData::get_vertex_array(MsqError &err) const 
  {
    if (vertexArray==0) 
      err.set_msg("\nWARNING: no coordinates array defined.\n");
    return vertexArray;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::set_vertex_coordinates"
  /*! PatchData::set_vertex_coordinates(Vector3D coords, size_t index, MsqError &err)

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
#define __FUNC__ "PatchData::get_element_array" 
  /*! PatchData::get_element_array(MsqError &err)

      \brief return the PatchData elements information as a C array of integers,
      i.e. a connectivity array.

      This can only be used if the storage mode has been set to RAW_ARRAYS.
      See set_storage_mode().
    */
  inline MsqMeshEntity* PatchData::get_element_array(MsqError &err) const
  {
    if (elementArray==0) 
      err.set_msg("\nWARNING: no element array defined.\n");
    return elementArray;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::get_vertex_to_elem_offset" 
  /*! PatchData::get_vertex_to_elem_offset(MsqError &err)

      \brief 
    */
  inline size_t* PatchData::get_vertex_to_elem_offset(MsqError &err) const
  {
    if (vertexToElemOffset==NULL) 
      err.set_msg("\nWARNING: no vertex to element data available.\n");
    return vertexToElemOffset;
  }

#undef __FUNC__
#define __FUNC__ "PatchData::get_vertex_to_elem_array" 
  /*! PatchData::get_vertex_to_elem_array(MsqError &err)

      \brief 
    */
  inline size_t* PatchData::get_vertex_to_elem_array(MsqError &err) const
  {
    if (elemsInVertex == NULL) 
      err.set_msg("\nWARNING: no vertex to element data available.\n");
    return elemsInVertex;
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
#define __FUNC__ "PatchData::create_coords_memento"
  /*! \fn PatchData::create_coords_memento(MsqError &err)
   This function instantiate PatchDataCoordsMemento object and returns a pointer to it.
   The PatchDataCoordsMemento contains the current state of the PatchData coordinates.
   It can be used to restore the same PatchData object to those coordinates.

   It is the responsibility of the caller to discard the PatchDataCoordsMemento
   when not needed any more.
  */
  inline PatchDataCoordsMemento* PatchData::create_coords_memento(MsqError &err)
  {
    PatchDataCoordsMemento* memento = new PatchDataCoordsMemento;
    memento->originator = this;
    if (numVertices)
      memento->coords = new Vector3D[numVertices];
    memento->numVertices = numVertices;
    
    // Copy the coordinates
    for (int i=0; i<numVertices; ++i)
    {
      memento->coords[i] = vertexArray[i];
    }
    return memento;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchData::set_to_coords_memento"
  /*! \fn PatchData::set_to_coords_memento(PatchDataCoordsMemento* memento, MsqError &err)
   This function restores a PatchData object coordinates to a previous state hold in
   a PatchDataCoordsMemento object (see create_coords_memento() ).

   The function checks whether the memento originates from this particular PatchData object.
   The function does not destroy the memento object: this is the caller responsibility.
  */
  inline void PatchData::set_to_coords_memento(
    PatchDataCoordsMemento* memento,
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
    for (int i=0; i<numVertices; ++i)
    {
      vertexArray[i] = memento->coords[i];
    }
  }
  
} // namespace


#endif
