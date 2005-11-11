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
#include "CornerTag.hpp"
#include "TargetMatrix.hpp"


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
  class Mesh;
  class TargetMatrix;
//   class SimplifiedGeometryEngine;
  
  /*! \class PatchData
    Contains all the mesh information necessary for
    one iteration of the optimization algorithms over a
    local mesh patch. */
  class PatchData
  {
  public:    
      // Constructor/Destructor
    MESQUITE_EXPORT PatchData();
    MESQUITE_EXPORT ~PatchData();
    
    
      /**\brief For use by testing code -- create patch explicitly
       * 
       * Create a patch containing elements of the same type and
       * without any higher-order nodes.
       *
       *\param num_vertex   Number of vertices in patch
       *\param vtx_coords   Array of vertex coords.  Length must be 3*num_vertex
       *\param type         Element type
       *\param connectivity Element connectivity, specified as a list
       *                    of vertex numbers, beginning with zero.
       *\param vertex_fixed_flags Optional array to specify which vertices
       *                    are to be marked as fixed.  If not specified,
       *                    no vertices are fixed.     
       */
    void fill( size_t num_vertex, const double* vtx_coords,
               size_t num_elem, EntityTopology type, 
               const size_t* connectivity,
               const bool* vertex_fixed_flags,
               MsqError& err );
    
      /**\brief For use by testing code -- create patch explicitly
       * 
       * Create a patch containing elements without any higher-order nodes.
       *
       *\param num_vertex   Number of vertices in patch
       *\param vtx_coords   Array of vertex coords.  Length must be 3*num_vertex
       *\param elem_types   The type of each element
       *\param connectivity Element connectivity, specified as a list
       *                    of vertex numbers, beginning with zero.
       *\param vertex_fixed_flags Optional array to specify which vertices
       *                    are to be marked as fixed.  If NULL,
       *                    no vertices are fixed.     
       */
    void fill( size_t num_vertex, const double* vtx_coords,
               size_t num_elem, const EntityTopology* elem_types,
               const size_t* connectivity,
               const bool* vertex_fixed_flags,
               MsqError& err );
    
      /**\brief For use by testing code -- create patch explicitly
       * 
       * Most general form of fill function.  Works for polygons, 
       * elements with higher-order nodes, etc.
       *
       *\param num_vertex   Number of vertices in patch
       *\param vtx_coords   Array of vertex coords.  Length must be 3*num_vertex
       *\param elem_types   The type of each element
       *\param vertex_per_elem The length of the connectivity list for each element.
       *\param connectivity Element connectivity, specified as a list
       *                    of vertex numbers, beginning with zero.
       *\param vertex_fixed_flags Optional array to specify which vertices
       *                    are to be marked as fixed.  If NULL,
       *                    no vertices are fixed.     
       */
    void fill( size_t num_vertex, const double* vtx_coords,
               size_t num_elem, const EntityTopology* elem_types,
               const size_t* vertex_per_elem,
               const size_t* elem_connectivity,
               const bool* vertex_fixed_flags,
               MsqError& err );
               
      /**\brief Create global patch
       *
       * Create a global patch - mesh should be initialized first.
       */
    void fill_global_patch( MsqError& err );
    
      /**\brief Iterate over mesh, creating element-on-vertex patches
       *
       * Get patch of elements adjacent to next vertex returned by
       * internal vertex iterator.  Mesh must be initialized before
       * calling this method.
       *\param num_layers Depth of adjacent elements to put in patch.
       *\param skip_fixed_vertices If true, vertices designated as fixed
       *         by the application or marked as culled (MSQ_SOFT_FIXED)
       *         will be skipped when iterating over all vertices.
       *\param free_vertex_index_out Output.  Index of the non-fixed
       *                   vertex (the 'current' vertex from the iterator)
       */
    bool get_next_vertex_element_patch( int num_layers,
                                        bool skip_fixed_vertices,
                                        size_t& free_vertex_index_out,
                                        MsqError& err );
    
      /**\brief Iterate over mesh, creating single-element patches
       *
       * Fill patch with next element returned by
       * internal element iterator.  Mesh must be initialized before
       * calling this method.
       */
    bool get_next_element_patch( MsqError& err );
   
     /**\brief Reset internal mesh iterator(s)
      *
      * Reset VeretxIterator used by \ref get_next_vertex_element_patch
      * and ElementIterator used by \ref get_next_element_patch
      */
    void reset_iterators();


  private:
      //! Doesn't allow PatchData to be copied implicitly.
      //! Mementos such as PatchDataVerticesMemento should be used when necessary. 
    PatchData(const PatchData &pd);
      //! Doesn't allow a PatchData object to be assigned to another.
      //! Mementos such as PatchDataVerticesMemento should be used when necessary. 
    PatchData& operator=(const PatchData &pd);
    
  public:

    enum ComputedInfo {
      MIN_UNSIGNED_AREA = 0, //!< minimum volume or area out of all elements in the patch
      MAX_UNSIGNED_AREA, //!< maximum volume or area out of all elements in the patch
      MIN_EDGE_LENGTH, //!< minimum edge length in the patch
      MAX_EDGE_LENGTH, //!< maximum edge length in the patch
      MINMAX_SIGNED_DET2D, //!< minimum and maximum corner area out of all elements in the patch
      MINMAX_SIGNED_DET3D, //!< minimum and maximum corner volume out of all elements in the patch
      AVERAGE_DET3D, //!< average corner determinant out of all elements in the patch
      MAX_COMPUTED_INFO_ENUM
    };

    //! This function clears the patch information such as maximum volume, etc ... 
    void clear_computed_info() { haveComputedInfos = 0; }
    
    bool have_computed_info( ComputedInfo info ) const
      { return 0 != (haveComputedInfos&(1<<info)); }
    
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

      //! Removes data
    void clear();
      //! Reorders the mesh data 
    void reorder();

      //! number of vertices in the patch. 
    size_t num_vertices() const
      { return numCornerVertices;}
      //! number of elements in the Patch.
    size_t num_elements() const
      { return elementArray.size(); }
      //! number of elements corners in the Patch. 
    size_t num_corners() ;
      /** Get number of nodes (vertex + higher-order nodes) */
    size_t num_nodes() const
      { return vertexArray.size(); }

      //! Returns the number of elements in the current patch who are
      //! free to move.  This is a costly function, since we have to check
      //! the flags of all vertices in the patch.
    int num_free_vertices(MsqError &err) const;
    unsigned num_free_nodes( MsqError& err ) const;
    
      //! Returns a pointer to the start of the vertex array.
    const MsqVertex* get_vertex_array( MsqError& err ) const;
    MsqVertex* get_vertex_array(MsqError &err);
    
      //! Returns a pointer to the start of the element array.
    const MsqMeshEntity* get_element_array( MsqError& err ) const;
    MsqMeshEntity* get_element_array(MsqError &err);
    
    size_t* get_connectivity_array( )
      { return &elemConnectivityArray[0]; }
      
    Mesh::ElementHandle* get_element_handles_array( )
      { return &elementHandlesArray[0]; }
    
    Mesh::VertexHandle* get_vertex_handles_array()
      { return &vertexHandlesArray[0]; }
    
      //! Returns the start of the vertex->element array.
      //! For each vertex in the patch, this array holds
      //! the number of elements the vertex is attached to,
      //! followed by the indices of those elements.
    //const size_t* get_vertex_to_elem_array(MsqError &err);
      //! Returns the start of the vertex->element offset
      //! array (v2e_o).  For vertex i, v2e_o[i] is the
      //! index into the vertex->element array (v2e) where
      //! vertex i's data begins.  So, v2e[v2e_o[i]] gives
      //! you the number of elements vertex i is attached
      //! to, and v2e[v2e_o[i]+1] gives you the index of
      //! the first element attached to vertex i.
    //const size_t* get_vertex_to_elem_offset(MsqError &err);
    
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
    
      /** Get the indices of elements adjacent to the specified vertex,
       *  and having the specified dimension */
    void get_vertex_element_indices( size_t vertex_index,
                                     unsigned element_dimension,
                                     msq_std::vector<size_t>& elem_indices,
                                     MsqError& err );
    
      /*! Get indices of elements attached to specified vertex */
    size_t* get_vertex_element_adjacencies( size_t vertex_index,
                                            size_t& array_len_out,
                                            MsqError& err );
    
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
    void generate_vertex_to_element_data();

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
    bool domain_set() const
    { return 0 != myDomain; }
    
      /*! Get the normal of the surface for a given vertex.
          Normal is returned in Vector3D &surf_norm.  If the normal cannot
          be determined, or if the underlying domain is not a surface,
          the normal will be set to (0,0,0).
          Check PatchData::domain_set() is not false first.
      */
    void get_domain_normal_at_vertex(size_t vertex_index, bool normalize,
                                     Vector3D &surf_norm,
                                     MsqError &err) ;
    
      /*! Get the normal to the domain at the centroid (projected to the
          domain) of a given element.
          Normal is returned in Vector3D &surf_norm.  If the normal cannot
          be determined, or if the underlying domain is not a surface,
          the normal will be set to (0,0,0).
          Check PatchData::domain_set() is not false first.
      */
    void get_domain_normal_at_element(size_t elem_index, Vector3D &surf_norm,
                                      MsqError &err) const;

      /** Get surface normal at a point where the surface is the
       *  domain of an element and the point is the location of 
       *  one of the element corners.
       */
    //void get_domain_normal_at_corner( size_t elem_index,
    //                                  size_t elem_corner,
    //                                  Vector3D& normal_out,
    //                                  MsqError& err ) const;
      /** Get surface normals at element corners.
       *  normals_out must be of sufficient size to hold
       *  the normals of all the corners.
       **/
    void get_domain_normals_at_corners( size_t element_index,
                                        Vector3D normals_out[],
                                        MsqError& err ) ;
                                        

      //! Alternative signature. Same functionality.
    void get_domain_normal_at_element(MsqMeshEntity* elem_ptr,
                                      Vector3D &surf_norm, MsqError &err) const 
    { get_domain_normal_at_element(size_t(elem_ptr-&(elementArray[0])), surf_norm, err); }
    
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
    PatchDataVerticesMemento* create_vertices_memento( MsqError &err,
                                                       bool include_higher_order = false );
    
      //! reinstantiates a memento to holds the current
      //! state of the PatchData coordinates. Improves memory management.
    void recreate_vertices_memento( PatchDataVerticesMemento* memento, 
                                    MsqError &err,
                                    bool include_higher_order = false );
    
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
        GLOBAL_PATCH,              /*!< Fills PatchData with all elements and vertices
                                    contained in all the meshes of the MeshSet. */
        ELEMENT_PATCH //!< fills PatchData with one element at a time.
      };

    //! Sets the originating meshSet. This is normally done in MeshSet::get_next_patch().
    //! This function is only for tests purposes. 
    void set_mesh(Mesh* ms);
    
    //! Returns the originating meshSet.
    Mesh* get_mesh() const
      { return myMesh; }
      
    void set_domain( MeshDomain* dm );
    
    MeshDomain* get_domain() const
      { return myDomain; }
    
    //! Target matrix data
    CornerTag<TargetMatrix> targetMatrices;
    void clear_tag_data();
    
    
    //! Display the coordinates and connectivity information
    friend msq_stdio::ostream& operator<<( msq_stdio::ostream&, const PatchData& );
   
   private:
   
    VertexIterator* vertex_iterator( MsqError& err );
    ElementIterator* element_iterator( MsqError& err );

      /** Call after filling vertex handle and connectivity arrays to
       * finish initializing the PatchData.  Reorders vertex handles array
       * such that all higher-order nodes are at end of array, updates
       * element connectivity array appropriately, initalizes numCornerVertices,
       * and per-element vertex and node counts.
       *
       * NOTE:  If the patch contains higher-order elements, this function
       *        will re-order the nodes in the vertex array. Do *NOT* assume
       *        vertex indices are the same after calling this function!
       *
       * NOTE:  This function expects the following data to be initalized:
       *         vertexHandlesArray
       *         elemConnectivityArray
       *         the topology type for all elements in elementArray
       *        The function assumes the following data has not been
       *        initialized and therefore does not need to be updated:
       *         vertexArray
       *
       * \param elem_offset_array Offset into connectivity array for each element
       */
    void initialize_data( size_t* elem_offset_array, MsqError& err );
   
      /** Code common to misc. methods for populating patch data.
       *  Calls \ref initialize_data, sets element types and
       *  sets vertex coordinates using vertexHandlesArray and 
       *  myMesh.
       *
       *  Note: The following data members must be initialized
       *  BEFORE CALLING THIS METHOD:
       *  - vertexHandlesArray
       *  - elemConnectivityArray
       *  - myMesh
       */
    void initialize_patch( EntityTopology* elem_type_array,
                           size_t* elem_offset_array,
                           MsqError& err );
    
      /** Code common to misc. methods for populating patch data.
       *  Remove duplicates from an array of handles.
       *\param handles    The array of handles to uniquify.
       *\param count      As input, the lenght of the \ref handles
       *                  array.  As output, the number of unique
       *                  handles remaining in the array.
       *\param index_map  If non-null, this must be an array of the
       *                  same length as the handles array.  If this
       *                  array is passed, the entry cooresponding
       *                  to each handle in the input \ref handles array will
       *                  be set to the index of that handle in the output
       *                  array.
       */
    static void make_handles_unique( Mesh::EntityHandle* handles,
                                     size_t& count,
                                     size_t* index_map = 0 );
    
      /*\brief Note that the passed info has been calculated and stored */
    void note_have_info( ComputedInfo info )
      { haveComputedInfos |= (1<<info); }
      
      /*\brief Update cached domain normal data */
    void update_cached_normals( MsqError& );

    Mesh* myMesh;              //!< The Mesh used to fill this PatchData [may be NULL]
    MeshDomain* myDomain;      //!< The geometric domain of the mesh [may be NULL]
    
    VertexIterator* vertexIterator;   //!< Current vertex in Mesh [may be NULL]
    ElementIterator* elementIterator; //!< Current element in Mesh [may be NULL]
    
      //! Cached data for vertices in \ref vertexHandlesArray,
      //! or vertex data for a temporary patch.
    msq_std::vector<MsqVertex> vertexArray;
      //! The list of handles for the vertices in this patch
      //! May be empty if \ref myMesh is NULL
    msq_std::vector<Mesh::VertexHandle> vertexHandlesArray;
      //! The number of vertices in \ref vertexArray that are
      //! corner vertices (not mid-nodes.)  This is the offset 
      //! in \ref vertexArray beginning at which the mid-side
      //! nodes are stored.
    size_t numCornerVertices;
      //! Cached data for elements in \ref elementHandlesArray
      //! or element data for a temporary patch.
    msq_std::vector<MsqMeshEntity> elementArray;
      //! The hist of handles for elements in this patch.
      //! May be empty if \ref myMesh is NULL
    msq_std::vector<Mesh::ElementHandle> elementHandlesArray;
      //! Element connectivity data.  The concatenation of the 
      //! connectivity list of each element in \ref elementArray.
      //! Each element in \ref elementArray has a pointer into
      //! this array at the correct offset for that element's connectivity data.
    msq_std::vector<size_t> elemConnectivityArray;
      //! The concatenation of the adjacency lists of all the vertices
      //! in \ref vertexArray.  Each value in the array is an index into 
      //! \ref elementArray indicating that the corresponding element uses
      //! the vertex.  May be empty if vertex adjacency data has not been
      //! requested.
    msq_std::vector<size_t> vertAdjacencyArray;
      //! This array is indexed by vertex indices and specifies the
      //! offset in \vertAdjacencyArray at which the adjacency list
      //! for the corresponding vertex begins.  May be empty if vertex 
      //! adjacency data has not been requested.
    msq_std::vector<size_t> vertAdjacencyOffsets;
      //! This array is indexed with vertex indices and contains a 
      //! pointer into \ref normalData at which the cached domain normal
      //! data for the vertex begins.  If the DOF in the domain is 2
      //! for the vertex, then the vertex has a single normal.  If
      //! the DOF is other than 2, the vertex has a normal for each
      //! attached two-dimensional element.  If normalData is not empty
      //! but this array is, then it is assumed that all vertices have a
      //! DOF of two and the single vertex normal is at an offset of the
      //! vertex index into the normalData array.
    msq_std::vector<Vector3D*> vertexNormalPointers;
      //! Storage space for cached domain normal data.  Pointers in
      //! \ref vertexNormalPointers point into this list.
    msq_std::vector<Vector3D> normalData;
      //! Storage space for cached domain DOF for vertices.  IF
      //! a domain exists and \ref normalData is not empty, but
      //! this array is, it may be assumed that all vertices have
      //! have a DOF == 2.
    msq_std::vector<unsigned short> vertexDomainDOF;
    
      // Arrays in which to store temporary data
      // (avoids reallocation of temp space)
    msq_std::vector<size_t> offsetArray;
    msq_std::vector<unsigned char> byteArray;
    
      // Patch Computed Information (maxs, mins, etc ... )
    double computedInfos[MAX_COMPUTED_INFO_ENUM];
      // Bit map indicating which values in \ref computedInfos
      // are valud (which values have been calculated.)
    unsigned haveComputedInfos;

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
    vertexArray.clear();
    vertexHandlesArray.clear();
    elementArray.clear();
    elementHandlesArray.clear();
    elemConnectivityArray.clear();
    vertAdjacencyArray.clear();
    vertAdjacencyOffsets.clear();
    vertexNormalPointers.clear();
    normalData.clear();
    vertexDomainDOF.clear();
    numCornerVertices = 0;
    haveComputedInfos = 0;
    myMesh = 0;
    myDomain = 0;
    delete vertexIterator;
    vertexIterator = 0;
    delete elementIterator;
    elementIterator = 0;
  }
  
  
  

  /*! \fn PatchData::get_vertex_array(MsqError &err) const 

  \brief Returns an array of all vertices in the PatchData.
  */
  inline const MsqVertex* PatchData::get_vertex_array(MsqError &err) const 
  {
    if (vertexArray.empty()) 
      MSQ_SETERR(err)( "No vertex array defined", MsqError::INVALID_STATE );
    return &vertexArray[0];
  }
  inline MsqVertex* PatchData::get_vertex_array(MsqError &err) 
  {
    if (vertexArray.empty()) 
      MSQ_SETERR(err)( "No vertex array defined", MsqError::INVALID_STATE );
    return &vertexArray[0];
  }
  
  /*! \fn PatchData::get_element_array(MsqError &err) const 

  \brief Returns the PatchData element array.
  */
  inline const MsqMeshEntity* PatchData::get_element_array(MsqError &err) const
  {
    if (elementArray.empty()) 
      MSQ_SETERR(err)( "No element array defined", MsqError::INVALID_STATE );
    return &elementArray[0];
  }
  inline MsqMeshEntity* PatchData::get_element_array(MsqError &err)
  {
    if (elementArray.empty()) 
      MSQ_SETERR(err)( "No element array defined", MsqError::INVALID_STATE );
    return &elementArray[0];
  }
  
  /*! \fn PatchData::set_vertex_coordinates(const Vector3D &coords, size_t index, MsqError &err)

  \brief set the coordinates of a vertex in the raw array
  */
  inline void PatchData::set_vertex_coordinates(const Vector3D &coords,
                                                size_t index,
                                                MsqError &err) 
  {
    if (index >= vertexArray.size()) {
      MSQ_SETERR(err)( "Index bigger than numVertices.", MsqError::INVALID_ARG );
      return;
    }
    
    vertexArray[index] = coords;
  }
  
  
  /*! \fn PatchData::get_vertex_to_elem_offset(MsqError &err) const 
   */  
  //inline const size_t* PatchData::get_vertex_to_elem_offset(MsqError &/*err*/)
  //{
  //    // Make sure we've got the data
  //  if (vertAdjacencyOffsets.empty())
  //  {
  //    generate_vertex_to_element_data();
  //  }
  //  return &vertAdjacencyOffsets[0];
  //}
  
  /*! \fn PatchData::get_vertex_to_elem_array(MsqError &err) const 
   */
  //inline const size_t* PatchData::get_vertex_to_elem_array(MsqError &/*err*/) 
  //{
  //    // Make sure we've got the data
  //  if (vertAdjacencyArray.empty())
  //  {
  //    generate_vertex_to_element_data();
  //  }
  //  return &vertAdjacencyArray[0];
  //}
  
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
    return vertex - &vertexArray[0];
  }
  
  inline size_t PatchData::get_element_index(MsqMeshEntity* element)
  {
    return element - &elementArray[0];
  }

  
  /*! \fn PatchData::create_vertices_memento(MsqError &err)
    This function instantiate PatchDataVerticesMemento object and returns a pointer to it.
    The PatchDataVerticesMemento contains the current state of the PatchData coordinates.
    It can be used to restore the same PatchData object to those coordinates.

    It is the responsibility of the caller to discard the PatchDataVerticesMemento
    when not needed any more.
  */
  inline PatchDataVerticesMemento* PatchData::create_vertices_memento(MsqError& /*err*/,
                                                                      bool include_higher_order)
  {
    size_t num_verts = include_higher_order ? num_nodes() : num_vertices();
    PatchDataVerticesMemento* memento = new PatchDataVerticesMemento;
    memento->originator = this;
    if (num_verts)
      memento->vertices = new MsqVertex[num_verts];
    memento->numVertices = num_verts;
    memento->arraySize = num_verts;
     
      // Copy the coordinates
    msq_stdc::memcpy(memento->vertices, &vertexArray[0], num_verts*sizeof(MsqVertex) );
    
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
                                                   MsqError& /*err*/,
                                                   bool include_higher_order)
  {
    size_t num_verts = include_higher_order ? num_nodes() : num_vertices();
    memento->originator = this;
    
    if ( num_verts > memento->arraySize
         || num_verts < memento->arraySize/10)
    {
      delete[] memento->vertices;
        // Allocate the new array
      memento->vertices = new MsqVertex[num_verts];
      memento->arraySize = num_verts;
    }
    
      // Copy the coordinates
    msq_stdc::memcpy(memento->vertices, &vertexArray[0],num_verts*sizeof(MsqVertex) );
    
    memento->numVertices = num_verts;
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
    
    if (memento->numVertices != num_vertices() &&
        memento->numVertices != num_nodes())
    {
      MSQ_SETERR(err)("Unable to restore patch coordinates.  Number of "
                      "vertices in PatchData has changed.",
                      MsqError::INVALID_STATE);
      return;
    }
    
      // copies the memento array into the PatchData array.
    msq_stdc::memcpy(&vertexArray[0], memento->vertices, memento->numVertices*sizeof(MsqVertex) );
  }
      
} // namespace


#endif
