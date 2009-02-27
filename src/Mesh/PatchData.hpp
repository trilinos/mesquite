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
#include "TargetMatrix.hpp"
#include "MsqError.hpp"
#include "Settings.hpp"

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
  class ExtraData;
  class PatchDataVerticesMemento;
  class Mesh;
  class Settings;
  class MappingFunction;
  class MappingFunction2D;
  class MappingFunction3D;
  
  /*!
    Contains all the mesh information necessary for
    one iteration of the optimization algorithms over a
    local mesh patch. */
  class PatchData 
  {
  public:    
      // Constructor/Destructor
    MESQUITE_EXPORT PatchData();
    MESQUITE_EXPORT ~PatchData();
    
    void attach_settings( const Settings* settings )
      { mSettings = settings; }
    const Settings* settings() const
      { return mSettings; }
    
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
    
    
    void set_mesh_entities( 
                   msq_std::vector<Mesh::ElementHandle>& patch_elems,
                   msq_std::vector<Mesh::VertexHandle>& free_vertices,
                   MsqError& err );
    
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
    
    void get_minmax_edge_length(double& min, double& max) const;
    
    //! Returns delta based on the minimum and maximum corner determinant over all elements in the patch
    //! This information is stored in the patch and should not decrease performance
    //! when used properly. See also PatchData::clear_computed_info() .
    double get_barrier_delta(MsqError &err); 

    //! Returns average corner determinant over all corners in the patch
    //! This information is stored in the patch and should not decrease performance
    //! when used properly. See also PatchData::clear_computed_info() .
//    double get_average_Lambda_3d(MsqError &err); 

      //! Removes data
    void clear();
      //! Reorders the mesh data 
    void reorder();

      //! number of vertices in the patch. 
    size_t num_nodes() const
      { return vertexArray.size();}
    size_t num_free_vertices() const
      { return numFreeVertices; }
    size_t num_slave_vertices() const
      { return numSlaveVertices; }
    size_t num_fixed_vertices() const
      { return num_nodes() - num_free_vertices() - num_slave_vertices(); }
      //! number of elements in the Patch.
    size_t num_elements() const
      { return elementArray.size(); }
    
    bool is_vertex_free(size_t index) const
      { return index < numFreeVertices; }
    bool is_vertex_not_free( size_t index ) const
      { return index >= numFreeVertices; }
    bool is_vertex_slave(size_t index) const
      { return index >= numFreeVertices && (index - numFreeVertices) < numSlaveVertices; }
    bool is_vertex_fixes(size_t index) const
      { return index >= numFreeVertices + numSlaveVertices; }
    
      //! number of element corners (number of vertex uses) in patch
    size_t num_corners() const;
    
      //! Returns a pointer to the start of the vertex array.
    const MsqVertex* get_vertex_array( MsqError& err ) const;
    //MsqVertex* get_vertex_array(MsqError &err);
    const MsqVertex* get_vertex_array() const { return &vertexArray[0]; }
    //MsqVertex* get_vertex_array()             { return &vertexArray[0]; }
    
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
    
    //MsqVertex& vertex_by_index(size_t index);
    const MsqVertex& vertex_by_index(size_t index) const;
    MsqMeshEntity& element_by_index(size_t index);
    const MsqMeshEntity& element_by_index(size_t index) const;
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

    void move_vertex( const Vector3D &delta, size_t index, MsqError &err);

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
    
      /*\brief Get domain normal at vertex location
       *
       * Get the normal of the domain associated with the passed
       * element handle at the location of the specified vertex.
       *\param vert_index  The index of the vertex in this PatchData
       *\param elem_handle The handle of the element passed to the domain 
       *\param normal_out  The resulting domain normal
       *\param err         Error flag.  Possible error conditions include:
       *                   invalid input data, no domain associated with
       *                   element, no domain at all, etc.
       */
    void get_domain_normal_at_vertex( size_t vert_index,
                                      Mesh::ElementHandle element,
                                      Vector3D& normal_out,
                                      MsqError& err );
    
      /*! Get the normal to the domain at the centroid (projected to the
          domain) of a given element.
          Normal is returned in Vector3D &surf_norm.  If the normal cannot
          be determined, or if the underlying domain is not a surface,
          the normal will be set to (0,0,0).
          Check PatchData::domain_set() is not false first.
      */
    void get_domain_normal_at_element(size_t elem_index, 
                                      Vector3D &surf_norm,
                                      MsqError &err);
    
      /** Get surface normals at element corners.
       *  normals_out must be of sufficient size to hold
       *  the normals of all the corners.
       **/
    void get_domain_normals_at_corners( size_t element_index,
                                        Vector3D normals_out[],
                                        MsqError& err ) ;
                                        
    void get_domain_normal_at_corner( size_t elemen_index,
                                      unsigned corner,
                                      Vector3D& normal,
                                      MsqError& err );
    
    void get_domain_normal_at_mid_edge( size_t element_index,
                                        unsigned edge_number,
                                        Vector3D& normal,
                                        MsqError& err );

      //! Alternative signature. Same functionality.
    void get_domain_normal_at_element(const MsqMeshEntity* elem_ptr,
                                      Vector3D &surf_norm, MsqError &err)  
    { get_domain_normal_at_element(size_t(elem_ptr-&(elementArray[0])), surf_norm, err); }
    
    void get_domain_normal_at_sample( size_t element_index,
                                      unsigned loc_dim,
                                      unsigned loc_num,
                                      Vector3D &surf_norm, MsqError &err)  
    {
      switch(loc_dim) {
        case 0:
          get_domain_normal_at_corner( element_index, loc_num, surf_norm, err );
          break;
        case 1:
          get_domain_normal_at_mid_edge( element_index, loc_num, surf_norm, err );
          break;
        case 2:
          assert(loc_num == 0);
          get_domain_normal_at_element( element_index, surf_norm, err );
          break;
        default:
          MSQ_SETERR(err)("Invalid dimension for surface element subentity.\n", MsqError::INVALID_ARG );
      }
    }

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
    
    //! Project gradient vector terms onto geometric domain
    void project_gradient( std::vector<Vector3D>& gradient, MsqError& err );

      //!Calculates the distance each vertex has moved from its original
      //!position as defined by the PatchDataVerticesMememnto.
    double get_max_vertex_movement_squared(PatchDataVerticesMemento* memento,
                                           MsqError &err);
    
      //! Updates the underlying mesh (the Mesquite::Mesh implementation) with
      //! new node coordinates and flag values.
    void update_mesh(MsqError &err);
    
      //! Calculate new location for all slave higher-order nodes using
      //! mapping function.  Called by update_mesh().
    void update_slave_node_coordinates( MsqError& err );
    void update_slave_node_coordinates( const size_t* elem_indices,
                                        size_t num_elem,
                                        MsqError& err );
    
      //!Remove the soft_fixed flag from all vertices in the patch.
    void set_all_vertices_soft_free(MsqError &err);
      //!Add a soft_fixed flag to all vertices in the patch.
    void set_all_vertices_soft_fixed(MsqError &err);
      //!Add a soft_fixed flag to all free vertices in the patch.
    void set_free_vertices_soft_fixed(MsqError &err);

      //!Mark vertex as culled (soft fixed)
    void set_vertex_culled( size_t vtx_index )
      { vertexArray[vtx_index].flags() |= MsqVertex::MSQ_CULLED; }
      //!Mark vertex as culled (soft fixed)
    void clear_vertex_culled( size_t vtx_index )
      { vertexArray[vtx_index].flags() &= ~MsqVertex::MSQ_CULLED; }
      //! check if vertex is culled
    int check_vertex_culled( size_t vtx_index ) const
      { return vertexArray[vtx_index].get_flags() | MsqVertex::MSQ_CULLED; }
    
      //! Fills a PatchData with the elements attached to a center vertex.
      //! Note that all entities in the sub-patch are copies of the entities
      //! in 'this' patch.  As such, moving a vertex in the sub-patch
      //! won't move the corresponding vertex in the source patch.  Also,
      //! calling 'update_mesh()' on the sub-patch WILL modify the TSTT
      //! mesh, but the source patch won't see the changes.
    void get_subpatch(size_t center_vertex_index,
                      unsigned num_adj_elem_layers,
                      PatchData &pd_to_fill,
                      MsqError &err);
    
    void get_free_vertex_coordinates( msq_std::vector<Vector3D>& coords_out ) const;
    
      //! Creates a memento that holds the current
      //! state of the PatchData coordinates. 
    PatchDataVerticesMemento* create_vertices_memento( MsqError &err );
    
      //! reinstantiates a memento to holds the current
      //! state of the PatchData coordinates. Improves memory management.
    void recreate_vertices_memento( PatchDataVerticesMemento* memento, 
                                    MsqError &err );
    
    //! Restore the PatchData coordinates to the state
    //! contained in the memento.
    void set_to_vertices_memento(PatchDataVerticesMemento* memento,
                                 MsqError &err);

    //! Sets the originating meshSet. This is normally done in MeshSet::get_next_patch().
    //! This function is only for tests purposes. 
    void set_mesh(Mesh* ms);
    
    //! Returns the originating meshSet.
    Mesh* get_mesh() const
      { return myMesh; }
      
    void set_domain( MeshDomain* dm );
    
    MeshDomain* get_domain() const
      { return myDomain; }
    
    
    const Settings* get_settings() const
      { return mSettings; }
    const MappingFunction* get_mapping_function( EntityTopology type ) const
      { return mSettings->get_mapping_function(type); }
    const MappingFunction2D* get_mapping_function_2D( EntityTopology type ) const
      { return mSettings->get_mapping_function_2D(type); }
    const MappingFunction3D* get_mapping_function_3D( EntityTopology type ) const
      { return mSettings->get_mapping_function_3D(type); }
    
    void get_sample_location( size_t element_index,
                              unsigned sample_dim,
                              unsigned sample_num,
                              Vector3D& result,
                              MsqError& err ) const;  
    
    //! This function produces an unsigned integer containing a
    //! bit for edge, face, and 1 for the region if the element 
    //! is a volume element.  The bits are ordered from lsb to msb
    //! in the same order as the corresponding mid-nodes would appear in
    //! the canonical ordering of an element with all such mid-nodes
    //! present.  (E.g. If the specified input element 
    //! has has HEXAHEDRON topology, then the result will have a
    //! bit for each of the 19 non-corner nodes in a Hex27, reglardless of
    //! how many nodes the actual indicated element has.)  Each bit
    //! will be 1 if the corresponding mid-node is a) present in
    //! this particular element and b) not a slave node.  Otherwise
    //! the bit will be 0.
    unsigned higher_order_node_bits( size_t elem_idx ) const;  
    
    //! Display the coordinates and connectivity information
    friend msq_stdio::ostream& operator<<( msq_stdio::ostream&, const PatchData& );
   
   private:
   
    /* allow access to the following to functions */
    friend class Mesquite::ExtraData;    
    /**\brief Attach an ExtraData object to this PatchData */
    bool attach_extra_data( ExtraData* data );
    /**\brief Remove an ExtraData object from this PatchData */
    bool remove_extra_data( ExtraData* data );

    /**\brief notify all attached ExtraData instances that patch contents have changed */
    void notify_new_patch( );
    /**\brief notify all attached ExtraData instances that a subpatch is being initalized 
     *\param sub_patch  The new, already populated subpatch
     *\param vertex_index_map For the i-th vertex in the subpatch, the
     *                  i-th entry in this list is the corresponding index in 
     *                  this patch.
     *\param element_index_map For the i-th element in the subpatch, the
     *                  i-th entry in this list is the corresponding index in 
     *                  this patch.
     */
    void notify_sub_patch( PatchData& sub_patch, 
                           const size_t* vertex_index_map,
                           const size_t* element_index_map,
                           MsqError& err );
    /**\brief notify all attached ExtraData instances that this patch is being destroyed */
    void notify_patch_destroyed();

      /** Call before initialize_data to change vertex_flags for 
       *  higher-order nodes to MSQ_DEPENDENT.
       */
    void enslave_higher_order_nodes( const size_t* element_offset_array,
                                     unsigned char* vertex_flags,
                                     MsqError& err ) const;

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
    void initialize_data( size_t* elem_offset_array, unsigned char* vertex_flags, MsqError& err );
    
      /** Code common to misc. methods for populating patch data.
       *  Remove duplicates from an array of handles.
       *\param handles    The array of handles to uniquify.
       *\param count      As input, the lenght of the #handles
       *                  array.  As output, the number of unique
       *                  handles remaining in the array.
       *\param index_map  If non-null, this must be an array of the
       *                  same length as the handles array.  If this
       *                  array is passed, the entry cooresponding
       *                  to each handle in the input #handles array will
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
    
      //! Cached data for vertices in PatchData::vertexHandlesArray,
      //! or vertex data for a temporary patch.
    msq_std::vector<MsqVertex> vertexArray;
      //! The list of handles for the vertices in this patch
      //! May be empty if PatchData::myMesh is NULL
    msq_std::vector<Mesh::VertexHandle> vertexHandlesArray;
      //! The number of vertices in PatchData::vertexArray that are
      //! free vertices.  The vertex array is sorted such that
      //! free vertices are first in the array.  This value
      //! is therefore the offset at which slave vertices begin
      //! in the array.
    size_t numFreeVertices;
      //! The number of slave vertices in vertexArray.  The
      //! vertices are ordered such that all fixed vertices occur
      //! after the slave vertices.  Therefore the offset at which
      //! the fixed vertices begin is numFreeVertices+numSlaveVertices
    size_t numSlaveVertices;
      //! Cached data for elements in PatchData::elementHandlesArray
      //! or element data for a temporary patch.
    msq_std::vector<MsqMeshEntity> elementArray;
      //! The hist of handles for elements in this patch.
      //! May be empty if PatchData::myMesh is NULL
    msq_std::vector<Mesh::ElementHandle> elementHandlesArray;
      //! Element connectivity data.  The concatenation of the 
      //! connectivity list of each element in PatchData::elementArray.
      //! Each element in PatchData::elementArray has a pointer into
      //! this array at the correct offset for that element's connectivity data.
    msq_std::vector<size_t> elemConnectivityArray;
      //! The concatenation of the adjacency lists of all the vertices
      //! in PatchData::vertexArray.  Each value in the array is an index into 
      //! PatchData::elementArray indicating that the corresponding element uses
      //! the vertex.  May be empty if vertex adjacency data has not been
      //! requested.
    msq_std::vector<size_t> vertAdjacencyArray;
      //! This array is indexed by vertex indices and specifies the
      //! offset in \vertAdjacencyArray at which the adjacency list
      //! for the corresponding vertex begins.  May be empty if vertex 
      //! adjacency data has not been requested.
    msq_std::vector<size_t> vertAdjacencyOffsets;
      //! Index into normalData at which the normal for the corresponding
      //! vertex index is located. Only vertices constrained to a single
      //! domain with a topological dimension of 2 have a unique domain
      //! normal.
    msq_std::vector<unsigned> vertexNormalIndices;
      //! Storage space for cached domain normal data.  Pointers in
      //! PatchData::vertexNormalPointers point into this list.
    msq_std::vector<Vector3D> normalData;
      //! Storage space for cached domain DOF for vertices.  IF
      //! a domain exists and PatchData::normalData is not empty, but
      //! this array is, it may be assumed that all vertices have
      //! have a DOF == 2.
    //msq_std::vector<unsigned short> vertexDomainDOF;
    
      // Arrays in which to store temporary data
      // (avoids reallocation of temp space)
    msq_std::vector<size_t> offsetArray;
    msq_std::vector<unsigned char> byteArray;
    msq_std::vector<bool> bitMap;
    bool* boolArray;
    size_t boolArraySize;
    bool* get_bool_array( size_t size ) {
      if (boolArraySize < size)
        boolArray = (bool*)realloc( boolArray, (boolArraySize = size) );
      return boolArray;
    }
    
      // Patch Computed Information (maxs, mins, etc ... )
    double computedInfos[MAX_COMPUTED_INFO_ENUM];
      // Bit map indicating which values in PatchData::computedInfos
      // are valud (which values have been calculated.)
    unsigned haveComputedInfos;
    
    ExtraData* dataList;
    
    const Settings* mSettings;

  };
  
  void print_patch_data( const PatchData& pd );
  
  
  /*! \brief Contains a copy of the coordinates of a PatchData.

    Use PatchDataVerticesMemento when you want to change the coordinates
    of a PatchData object but also have the option to restore them.
    This class can only be instantiated through PatchData::create_vertices_memento().
  */
  class PatchDataVerticesMemento
  {
  public:
    void clear()
    {
      originator = 0;
      vertices.clear();
      normalData.clear();
    }
  private:
    // Constructor accessible only to originator (i.e. PatchData)
    friend class PatchData;
    PatchDataVerticesMemento() 
      : originator(0)
    {}
    
    PatchData* originator; //!< PatchData whose state is kept
    msq_std::vector<MsqVertex> vertices;
    msq_std::vector<Vector3D> normalData;
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
    vertexNormalIndices.clear();
    normalData.clear();
    //vertexDomainDOF.clear();
    numFreeVertices = 0;
    numSlaveVertices = 0;
    haveComputedInfos = 0;
    myMesh = 0;
    myDomain = 0;
  }
  
  
  

  /*! \brief Returns an array of all vertices in the PatchData.
  */
  inline const MsqVertex* PatchData::get_vertex_array(MsqError &err) const 
  {
    if (vertexArray.empty()) 
      MSQ_SETERR(err)( "No vertex array defined", MsqError::INVALID_STATE );
    return &vertexArray[0];
  }
  //inline MsqVertex* PatchData::get_vertex_array(MsqError &err) 
  //{
  //  if (vertexArray.empty()) 
  //    MSQ_SETERR(err)( "No vertex array defined", MsqError::INVALID_STATE );
  //  return &vertexArray[0];
  //}
  
  /*! \brief Returns the PatchData element array.
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
 
  /*! \brief set the coordinates of a vertex in the raw array
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
    
    if (numSlaveVertices) {
      size_t num_elem, *indices;
      indices = get_vertex_element_adjacencies( index, num_elem, err ); MSQ_ERRRTN(err);
      update_slave_node_coordinates( indices, num_elem, err ); MSQ_ERRRTN(err);
    }
  }
  inline void PatchData::move_vertex( const Vector3D &delta,
                                      size_t index,
                                      MsqError &err) 
  {
    if (index >= vertexArray.size()) {
      MSQ_SETERR(err)( "Index bigger than numVertices.", MsqError::INVALID_ARG );
      return;
    }
    
    vertexArray[index] += delta;
    
    if (numSlaveVertices) {
      size_t num_elem, *indices;
      indices = get_vertex_element_adjacencies( index, num_elem, err ); MSQ_ERRRTN(err);
      update_slave_node_coordinates( indices, num_elem, err ); MSQ_ERRRTN(err);
    }
  }
  
  //inline MsqVertex& PatchData::vertex_by_index(size_t index)
  //{
  //  return vertexArray[index];
  //}
  
  inline const MsqVertex& PatchData::vertex_by_index( size_t index ) const
    { return vertexArray[index]; }
  
  inline MsqMeshEntity& PatchData::element_by_index(size_t index)
  {
    return elementArray[index];
  }
  
  inline const MsqMeshEntity& PatchData::element_by_index( size_t index ) const
    { return elementArray[index]; }
  
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

  inline void PatchData::get_free_vertex_coordinates( msq_std::vector<Vector3D>& coords_out ) const
  {
    coords_out.resize( num_free_vertices() );
    msq_std::copy( vertexArray.begin(), vertexArray.begin()+num_free_vertices(), 
                   coords_out.begin() );
  }
    

  
  /*! 
    This function instantiate PatchDataVerticesMemento object and returns a pointer to it.
    The PatchDataVerticesMemento contains the current state of the PatchData coordinates.
    It can be used to restore the same PatchData object to those coordinates.

    It is the responsibility of the caller to discard the PatchDataVerticesMemento
    when not needed any more.
  */
  inline PatchDataVerticesMemento* PatchData::create_vertices_memento(MsqError& err)
  {
    PatchDataVerticesMemento* memento = new PatchDataVerticesMemento;
    recreate_vertices_memento( memento, err );
    if (MSQ_CHKERR(err)) {
      delete memento;
      return 0;
    }
    return memento;
  }
  
  /*! 
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
    
    size_t num_vtx = num_free_vertices() + num_slave_vertices();
    
    memento->vertices.resize( num_vtx );
    msq_std::copy( vertexArray.begin(), vertexArray.begin()+num_vtx, memento->vertices.begin() );
    
    int num_normal;
    if (normalData.empty())
      num_normal = 0;
    else if (vertexNormalIndices.empty()) 
      num_normal = num_vtx;
    else {
      num_normal = num_vtx;
      while (num_normal != 0 && vertexNormalIndices[--num_normal] >= normalData.size());
      if (num_normal == 0) {
        if (vertexNormalIndices[0] < normalData.size())
          num_normal = vertexNormalIndices[0] + 1;
      }
      else 
        num_normal = vertexNormalIndices[num_normal] + 1;
    }
    
    memento->normalData.resize( num_normal );
    msq_std::copy( normalData.begin(), normalData.begin()+num_normal, memento->normalData.begin() );
  }
  
  /*! 
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
    
    if (memento->vertices.size() != num_free_vertices()+num_slave_vertices())
    {
      MSQ_SETERR(err)("Unable to restore patch coordinates.  Number of "
                      "vertices in PatchData has changed.",
                      MsqError::INVALID_STATE);
      return;
    }
    
      // copies the memento array into the PatchData array.
    msq_std::copy( memento->vertices.begin(), memento->vertices.end(), vertexArray.begin() );
    msq_std::copy( memento->normalData.begin(), memento->normalData.end(), normalData.begin() );
  }
      
} // namespace


#endif
