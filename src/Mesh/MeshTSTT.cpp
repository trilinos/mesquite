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
/*!
  \file   MeshTSTT.cpp
  \brief  


  \author Jason Kraftcheck
  \date   2004-10-26
*/

#include "MeshTSTT.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "TSTTM.hh"
#include "TSTTM_EntityTopology.hh"
#include "sidl_cxx.hh"

namespace Mesquite
{

const char* const VERTEX_BYTE_TAG_NAME  = "MesquiteVertexByte";
const char* const VERTEX_FIXED_TAG_NAME = "MesquiteVertexFixed";

static msq_std::string process_tstt_error( TSTT::Error &tstt_err )
{
  msq_std::string str;
  msq_std::string result("TSTT ERROR: ");
  result += tstt_err.getNote();
  MSQ_DBGOUT(1) << "TSTT Error:" << msq_std::endl;
  MSQ_DBGOUT(1) << tstt_err.getNote() << msq_std::endl;
  tstt_err.getDescription(str);
  MSQ_DBGOUT(1) << str << msq_std::endl;
  MSQ_DBGOUT(1) << tstt_err.getTrace() << msq_std::endl;
  return result;
}

template <class T> static sidl::array<T> alloc_sidl_vector( size_t size )
{
  int32_t lower = 0;
  int32_t upper = size - 1;
  return sidl::array<T>::createRow( 1, &lower, &upper );
}

template <class T> static sidl::array<T> alloc_sidl_vector( size_t size, T init )
{
  sidl::array<T> result = alloc_sidl_vector(size);
  for (int32_t i = 0; i < (int32_t)size; ++i)
    result.set( i, init );
  return result;
}

template <class S, class T> static void copy_from_sidl( sidl::array<S>& source,
                                                        T* target )
{
  typename sidl::array<S>::iterator i = source.begin();
  for (; i != source.end(); ++i, ++target)
    *target = (T)*i;
}

/*************************************************************************
 *                          Iterator Definition
 ************************************************************************/

/** \brief Wrapper around TSTT interator interface */
class TSTTIterator : public EntityIterator
{
  private:
    TSTTM::Entity tsttMesh;        /**< TSTT mesh interface */
    void* tsttIter;                /**< TSTT iterator handle */
    bool notAtEnd;                 /**< Flag to mark if end of iterator had been reached */
    void* entityHandle;            /**< Current TSTT entity handle */

      /** 
       * Advance the iterator to the next entity.  Updates \ref entityHandle and 
       * \ref notAtEnd
       */
    inline void get_next_entity()
      { notAtEnd = tsttMesh.getNextEntIter( tsttIter, entityHandle ); }

  public:

      /**
       *\param mesh    The TSTT mesh interface instance
       *\param meshset The meshset to iterate over
       *\param type    TSTT type in \ref meshset to iterate over
       *\param topo    TSTT topology in \ref meshset to iterate over
       */
    TSTTIterator( TSTTM::Entity& mesh, 
                  void* meshset, 
                  TSTTM::EntityType type,
                  TSTTM::EntityTopology topo,
                  MsqError& err );
      /**\brief reset iterator */
    virtual void restart();
      /**\brief get current entity handle */
    virtual Mesh::EntityHandle operator*() const;
      /**\brief check if any remaining entity handles */
    virtual bool is_at_end() const;
      /**\biref step */
    virtual void operator++();
    
    virtual ~TSTTIterator();
};



TSTTIterator::TSTTIterator( TSTTM::Entity& mesh, void* meshset, 
                            TSTTM::EntityType type, TSTTM::EntityTopology topo,
                            MsqError& err )
  : tsttMesh(mesh), tsttIter(0), notAtEnd(true), entityHandle(0)
{
  try {
    tsttMesh.initEntIter( meshset, type, topo, tsttIter );
    get_next_entity();
  } 
  catch (TSTT::Error& tstt_err) {
    MSQ_SETERR(err)( process_tstt_error( tstt_err ), MsqError::INTERNAL_ERROR );
  }
}

TSTTIterator::~TSTTIterator()
{
  if (tsttIter)
    tsttMesh.endEntIter( tsttIter );
}

void TSTTIterator::restart()
{
  try {
    tsttMesh.resetEntIter( tsttIter );
    get_next_entity();
  }
  catch (TSTT::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}

Mesh::EntityHandle TSTTIterator::operator*() const
{
  return entityHandle;
}

bool TSTTIterator::is_at_end() const
{
  return !notAtEnd;
}

void TSTTIterator::operator++() 
{
  try {
    get_next_entity();
  }
  catch (TSTT::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}
       

/** \brief Iterate over a sidl array of TSTT entity handles */
class SIDLIterator : public EntityIterator
{
  private:
    sidl::array<void*>::const_iterator iter;
    const sidl::array<void*> array;
  
  public:

      /**\param array Array to iterate over */
    SIDLIterator( const sidl::array<void*>& a )
      : iter( a.begin() ), array( a )            {}

      /**\brief reset iterator */
    virtual void restart()                       { iter = array.begin(); }
    
      /**\brief get current entity handle */
    virtual Mesh::EntityHandle operator*() const { return *iter; }
    
      /**\brief check if any remaining entity handles */
    virtual bool is_at_end() const               { return iter == array.end(); }
    
      /**\biref step */
    virtual void operator++()                    { ++iter; }
};




/*************************************************************************
 *                            Mesh Interface 
 ************************************************************************/

  /** \brief Implementation of MeshTSTT
   */
  class MeshTSTTImpl : public MeshTSTT
  {
  public:

    MeshTSTTImpl(TSTTM::Mesh& tstt_mesh, MsqError &err);
    virtual ~MeshTSTTImpl();
    
      /** \brief set mesh to be smoothed.
       *
       * Set the mesh which Mesquite is to smooth.  Optionally
       * specify fixed vertices.
       * NOTE: If an active set is not specified, the default
       *       is to use the global set (the ENTIRE mesh.)
       *
       *\param element_set TSTT entity set handle for set containing
       *                  mesh elements for which quality is to be improved.
       *                  Any non-element entities in the set will be
       *                  ignored.
       *\param vertex_set TSTT entity set handle containing list of 
       *                  vertices which are not to be moved.  If NULL,
       *                  it is assumed that all vertices may be moved.
       */
    virtual void set_active_set( void* element_set, void* vertex_set, MsqError& );
    

      /**\brief Get dimension of vertex coordinates (2D vs. 3D). */
    virtual int get_geometric_dimension(Mesquite::MsqError &/*err*/);
    
    /** \brief get sizes for calling \ref get_all_mesh
     *
     * Get counts of entities in mesh.
     *
     *\param vertex_count  - Number of vertices connected to active mesh
     *\param element_count - Number of elements in active mesh
     *\param vertex_use_count - Number of vertex uses (sum of the length
     *                          of the connectivity list for all elements
     *                          in active.)
     */
    virtual void get_all_sizes( size_t& vertex_count,
                                size_t& element_count,
                                size_t& vertex_use_count,
                                MsqError& err );
    
    /** \brief Get entities and connectivity 
     *
     * Get vertex handles, element handles, and connectivty
     * for active mesh.  Use \ref get_all_sizes to determine
     * required array sizes.
     *
     *\param vert_array        Array to store vertex handles in
     *\param vert_len          Length of \ref vert_array
     *\param elem_array        Array to store element handles in
     *\param elem_len          Length of \ref elem_array
     *\param elem_conn_offsets Offsets into \ref elem_conn_indices at
     *                         which the connectivity data for each
     *                         element begins.  
     *\param offset_len        Length of \ref elem_conn_offsets.  Should
     *                         be \ref elem_len + 1.
     *\param elem_conn_indices Indices into \ref vert_array
     *\param index_len         Length of \ref elem_conn_indices.
     */
    virtual void get_all_mesh( VertexHandle*  vert_array, size_t vert_len,
                               ElementHandle* elem_array, size_t elem_len,
                               size_t* elem_conn_offsets, size_t offset_len,
                               size_t* elem_conn_indices, size_t index_len,
                               MsqError& err );
    
      /**\brief Create iterator for vertices in active set */
    virtual VertexIterator* vertex_iterator(MsqError &err);
    
      /**\brief Create iterator for elements in active set */
    virtual ElementIterator* element_iterator(MsqError &err);

      /**\brief Query "fixed" flag for a vertex */
    virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err);

      /**\brief Query "boundary" flag for an array of vertices */
    virtual void vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
                                 size_t num_vtx, MsqError &err);
    
      /**\brief Get vertex coordinates */
    virtual void vertices_get_coordinates(VertexHandle vert_array[],
                                  MsqVertex* const &coordinates,
                                  const size_t &num_vtx, MsqError &err);
      /**\brief Set vertex coordinates */
    virtual void vertex_set_coordinates(VertexHandle vertex,
                                 const Vector3D &coordinates, MsqError &err);
    
      /**\brief Set vertex mark */
    virtual void vertex_set_byte (VertexHandle vertex,
                            unsigned char byte, MsqError &err);
      /**\brief Set vertex mark */
    virtual void vertices_set_byte (VertexHandle *vert_array,
                              unsigned char *byte_array,
                              size_t array_size, MsqError &err);
    
      /**\brief Get vertex mark */
    virtual void vertex_get_byte(VertexHandle vertex,
                                 unsigned char *byte, MsqError &err);
      /**\brief Get vertex mark */
    virtual void vertices_get_byte(VertexHandle *vert_array,
                                   unsigned char *byte_array,
                                   size_t array_size, MsqError &err);
    
      /**\brief Get vertex adjacencies */
    virtual size_t vertex_get_attached_element_count(VertexHandle vertex, MsqError &err);
    
      /**\brief Get vertex adjacencies */
    virtual void vertex_get_attached_elements(VertexHandle vertex,
                                              ElementHandle* elem_array,
                                              size_t sizeof_elem_array,
                                              MsqError &err);
    
      /**\biref Get length of connectivity list */
    virtual size_t element_get_attached_vertex_count(ElementHandle elem,
                                                  MsqError &err);
    
/**\brief Get element connectivity in overly-complex CSR rep.
 *
 * Returns the vertices that are part of the topological definition of each
 * element in the "elem_handles" array.  When this function is called, the
 * following must be true:
 *   a) "elem_handles" points at an array of "num_elems" element handles.
 *   b) "vert_handles" points at an array of size "sizeof_vert_handles"
 *   c) "csr_data" points at an array of size "sizeof_csr_data"
 *   d) "csr_offsets" points at an array of size "num_elems+1"
 *      
 * When this function returns, adjacency information will be stored
 * in csr format:
 *    a) "vert_handles" stores handles to all vertices found in one
 *       or more of the elements.  Each vertex appears only
 *       once in "vert_handles", even if it is in multiple elements.
 *    b) "sizeof_vert_handles" is set to the number of vertex
 *       handles placed into "vert_handles".
 *    c) "sizeof_csr_data" is set to the total number of vertex uses (for
 *       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
 *       the two triangles share some vertices).
 *    c) "csr_offsets" is filled such that csr_offset[i] indicates the location
 *       of entity i's first adjacency in "csr_data".  The number of vertices
 *       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
 *       reason, csr_offsets[num_elems] is set to the new value of
 *       "sizeof_csr_data".
 *    d) "csr_data" stores integer offsets which give the location of
 *       each adjacency in the "vert_handles" array.
 *
 * As an example of how to use this data, you can get the handle of the first
 * vertex in element #3 like this:
 *   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ]
 *
 * and the second vertex of element #3 like this:
 *   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ]
 */ 
    virtual void elements_get_attached_vertices(ElementHandle *elem_handles,
                                                size_t num_elems,
                                                VertexHandle *vert_handles,
                                                size_t &sizeof_vert_handles,
                                                size_t *csr_data,
                                                size_t &sizeof_csr_data,
                                                size_t *csr_offsets,
                                                MsqError &err);
    
  
      /**\brief Return topology type enum for specified element */
    virtual EntityTopology element_get_topology(ElementHandle entity_handle,
                                           MsqError &err);
      /**\brief Return topology type enum for an array of elements */
    virtual void elements_get_topologies(ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err);
    
//**************** Memory Management ****************
      /**\brief no-op */ 
    virtual void release_entity_handles(EntityHandle *handle_array,
                                        size_t num_handles, MsqError &err);
    
      // Instead of deleting a Mesh when you think you are done,
      // call release().  In simple cases, the implementation could
      // just call the destructor.  More sophisticated implementations
      // may want to keep the Mesh object to live longer than Mesquite
      // is using it.
    virtual void release();

//*************** Tags  ***********

      /** \brief Create a tag
       *
       * Create a user-defined data type that can be attached
       * to any element or vertex in the mesh.  For an opaque or
       * undefined type, use type=BYTE and length=sizeof(..).
       *
       * \param tag_name  A unique name for the data object
       * \param type      The type of the data
       * \param length    Number of values per entity (1->scalar, >1 ->vector)
       * \param default_value Default value to assign to all entities - may be NULL
       * \return - Handle for tag definition 
       */
    virtual TagHandle tag_create( const msq_std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err);
     
      /** \brief Remove a tag and all corresponding data
       *
       * Delete a tag.
       */
    virtual void tag_delete( TagHandle handle, MsqError& err );
    
    
      /** \brief Get handle for existing tag, by name. */
    virtual TagHandle tag_get( const msq_std::string& name, 
                               MsqError& err );
     
      /** \brief Get properites of tag
       *
       * Get data type and number of values per entity for tag.
       * \param handle     Tag to get properties of.
       * \param name_out   Passed back tag name.
       * \param type_out   Passed back tag type.
       * \param length_out Passed back number of values per entity.
       */
    virtual void tag_properties( TagHandle handle,
                                 msq_std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err );
    
      /** \brief Set tag values on elements
       * 
       * Set the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err );

      /** \brief Set tag values on vertices
       * 
       * Set the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of node_array
       * \param node_array Array of vertices for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err );
    
    
      /** \brief Get tag values on elements
       * 
       * Get the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err );
    
      /** \brief Get tag values on vertices.
       * 
       * Get the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of vertices for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err );
  private:
      /** \brief Set tag values */
    void tag_set_data ( TagHandle handle,
                        size_t num_elems,
                        const EntityHandle* handle_array,
                        const void* tag_data,
                        MsqError& err );
    
      /** \brief Get tag values */
    void tag_get_data( TagHandle handle,
                       size_t num_elems,
                       const EntityHandle* handle_array,
                       void* tag_data,
                       MsqError& err );
    
    /** TSTT basic mesh interface instance */
    TSTTM::Mesh meshIFace;
    /** TSTT interface for per-entity queries */
    TSTTM::EntTag entIFace;
    /** TSTT interface for multi-entity (array) queries */
    TSTTM::ArrTag arrIFace;
    /** TSTT interface for modifying mesh */
    TSTTM::Modify modIFace;
    /** TSTT interface for entity set operations */
    TSTTM::EntSet setIFace;
    
    /** TSTTM entity set handle for fixed vertices */
    void* fixedVertices;
    /** TSTTM entity set handle for elements to improve */
    void* activeSet;
    
    /** Handle for tag used to hold vertex byte */
    TagHandle byteTag; 
    /** Tag was created in constructor */
    bool createdByteTag;
    /** Handle for tag used to hold vertex-fixed flag */
    TagHandle fixedTag;
    /** Fixed tag was created in constructor */
    bool createdFixedTag;
    
    /** Map TSTTM::EntityTopology to Mesquite::EntityTopology */
    EntityTopology topologyMap[TSTTM::EntityTopology_ALL_TOPOLOGIES];
 };

/*************************************************************************
 *                          Mesh Definition
 ************************************************************************/

MeshTSTT* MeshTSTT::create( TSTTM::Mesh& mesh, void* meshset, MsqError& err )
{
  MeshTSTT* result = new MeshTSTTImpl( mesh, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    result = 0;
  }
  result->set_active_set( meshset, 0, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    result = 0;
  }
  return result;
}

MeshTSTT::~MeshTSTT() {}


MeshTSTTImpl::MeshTSTTImpl(TSTTM::Mesh& tstt_mesh, Mesquite::MsqError& err) 
  : meshIFace(tstt_mesh), 
    fixedVertices(0), activeSet(0), 
    byteTag(0), createdByteTag(false),
    fixedTag(0), createdFixedTag(false)
{
    // Initialize topology map 
  
  const size_t mapsize = sizeof(topologyMap) / sizeof(Mesquite::EntityTopology);
  if (mapsize < TSTTM::EntityTopology_ALL_TOPOLOGIES)
  {
    MSQ_SETERR(err)("MeshTSTT needs to be updated for new TSTT element topologies.",
                    MsqError::INTERNAL_ERROR);
    return;
  }
  
  for (size_t i = 0; i <= TSTTM::EntityTopology_ALL_TOPOLOGIES; ++i)
    topologyMap[i] = Mesquite::MIXED;
  
  topologyMap[TSTTM::EntityTopology_TRIANGLE     ] = Mesquite::TRIANGLE;
  topologyMap[TSTTM::EntityTopology_QUADRILATERAL] = Mesquite::QUADRILATERAL;
  topologyMap[TSTTM::EntityTopology_TETRAHEDRON  ] = Mesquite::TETRAHEDRON;
  topologyMap[TSTTM::EntityTopology_HEXAHEDRON   ] = Mesquite::HEXAHEDRON;
  topologyMap[TSTTM::EntityTopology_PRISM        ] = Mesquite::PRISM;
  topologyMap[TSTTM::EntityTopology_PYRAMID      ] = Mesquite::PYRAMID;
  
  try {
      // Attempt to cast to single-entity query interface
    entIFace  = tstt_mesh;
    if (!entIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::EntTag",
                       MsqError::INVALID_STATE );
      return;
    } 
      // Attempt cast to multi-entity query interface
    arrIFace  = tstt_mesh;
    if (!arrIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::ArrTag",
                       MsqError::INVALID_STATE );
      return;
    }
      // Attempt cast to modify interface
    modIFace  = tstt_mesh;
    if (!arrIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::Modify",
                       MsqError::INVALID_STATE );
      return;
    }
      // Attempt cast to entity set interface
    setIFace  = tstt_mesh;
    if (!arrIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::EntSet",
                       MsqError::INVALID_STATE );
      return;
    }
    
      // Get/create tag for fixed flag
    try {
      bool def = false;
      entIFace.createTag( VERTEX_FIXED_TAG_NAME, 1, TSTT::TagValueType_BOOLEAN, &def, fixedTag );
      if (fixedTag)
        createdFixedTag = true;
    } 
    catch (TSTT::Error& tstt_err) {
      fixedTag = entIFace.getTagHandle( VERTEX_FIXED_TAG_NAME );
      if (!fixedTag)
        throw;
    }  
      // Double-check types incase tag already existed
    if (entIFace.getTagSize(fixedTag) != 1 || 
        entIFace.getTagSize(fixedTag) != TSTT::TagValueType_BOOLEAN)
    {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                       "Tag \"%s\" exists with invalid type/size", 
                       VERTEX_FIXED_TAG_NAME );
      return;
    }
    
      // Get/create tag for vertex byte
    try {
      char def = '\0';
      entIFace.createTag( VERTEX_BYTE_TAG_NAME, 1, TSTT::TagValueType_OPAQUE, &def, byteTag );
      if (byteTag)
        createdByteTag = true;
    } 
    catch (TSTT::Error& tstt_err) {
      fixedTag = entIFace.getTagHandle( VERTEX_FIXED_TAG_NAME );
      if (!fixedTag)
        throw;
    }  
      // Double-check types incase tag already existed
    if (entIFace.getTagSize(fixedTag) != 1 || 
        entIFace.getTagType(fixedTag) != TSTT::TagValueType_OPAQUE)
    {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                       "Tag \"%s\" exists with invalid type/size", 
                       VERTEX_FIXED_TAG_NAME );
      return;
    }

  }
  catch (TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    try {
      if (createdFixedTag)
        entIFace.destroyTag( fixedTag, false );
      if (createdByteTag)
        entIFace.destroyTag( byteTag, false );
    } 
    catch (...) {}
  }
}

Mesquite::MeshTSTTImpl::~MeshTSTTImpl() 
{
  try {
    if (createdFixedTag)
      entIFace.destroyTag( fixedTag, false );
    if (createdByteTag)
      entIFace.destroyTag( byteTag, false );
  } 
  catch (TSTT::Error& tstt_err ) {
    process_tstt_error(tstt_err);
    throw;
  }
}


void MeshTSTTImpl::set_active_set( void* elem_set, void* fixed_verts, MsqError& err )
{
  try {
    void* iter;
    void* handle;
  
      // Make sure elem_set contains only faces and regions
    entIFace.initEntIter( elem_set, TSTTM::EntityType_ALL_TYPES, TSTTM::EntityTopology_ALL_TOPOLOGIES, iter );
    while (entIFace.getNextEntIter( iter, handle ))
    {
      TSTTM::EntityType type = entIFace.getEntType( handle );
      if (type != TSTTM::EntityType_FACE && type != TSTTM::EntityType_REGION)
      {
        MSQ_SETERR(err)("Active set must contain only FACEs and/or REGIONs",
                        MsqError::INVALID_ARG);
        return;
      }
    }
  
  
      // Clear old fixed flags
    if (fixedVertices)
    {
      entIFace.initEntIter( fixedVertices, TSTTM::EntityType_VERTEX, TSTTM::EntityTopology_POINT, iter );
      while (entIFace.getNextEntIter( iter, handle ))
        entIFace.setBoolData( handle, fixedTag, false );
      entIFace.endEntIter( iter );
    }
    
    fixedVertices = fixed_verts;
    
      // Set new fixed flags
    if (fixedVertices)
    {
      entIFace.initEntIter( fixedVertices, TSTTM::EntityType_VERTEX, TSTTM::EntityTopology_POINT, iter );
      while (entIFace.getNextEntIter( iter, handle ))
        entIFace.setBoolData( handle, fixedTag, false );
      entIFace.endEntIter( iter );
    }
    
    activeSet = elem_set;
  }
  catch (TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}    
      
  

// Returns whether this mesh lies in a 2D or 3D coordinate system.
int MeshTSTTImpl::get_geometric_dimension(Mesquite::MsqError &err)
{
  try {
    return meshIFace.getGeometricDim();
  }
  catch (TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}
    
    
// Returns a pointer to an iterator that iterates over the
// set of all vertices in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If vertices are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
VertexIterator* MeshTSTTImpl::vertex_iterator(MsqError& err)
{
  try {
  
    sidl::array<void*> handles;
    sidl::array<int> offsets, indices;
    int num_handles, num_offsets, num_indices;
    meshIFace.getAdjEntities( activeSet, 
                              TSTTM::EntityType_ALL_TYPES,
                              TSTTM::EntityTopology_ALL_TOPOLOGIES,
                              TSTTM::EntityType_VERTEX,
                              handles, num_handles,
                              offsets, num_offsets,
                              indices, num_indices );
                              
    handles = alloc_sidl_vector<void*>(num_handles);
    offsets = alloc_sidl_vector<int>(num_offsets);
    indices = alloc_sidl_vector<int>(num_indices);
    meshIFace.getAdjEntities( activeSet, 
                              TSTTM::EntityType_ALL_TYPES,
                              TSTTM::EntityTopology_ALL_TOPOLOGIES,
                              TSTTM::EntityType_VERTEX,
                              handles, num_handles,
                              offsets, num_offsets,
                              indices, num_indices );
    
    return new SIDLIterator( handles );
    
  }
  catch (TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}
    
// Returns a pointer to an iterator that iterates over the
// set of all top-level elements in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If elements are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
ElementIterator* MeshTSTTImpl::element_iterator(MsqError &err)
{
  // Use array instead of TSTT iterator for efficiency
  //
  //ElementIterator* result = new TSTTIterator( entIFace, 
  //                                            activeSet,
  //                                            TSTTM::EntityType_ALL_TYPES,
  //                                            TSTTM::EntityTopology_ALL_TOPOLOGIES,
  //                                            err );
  //if (MSQ_CHKERR(err))
  //{
  //  delete result;
  //  return 0;
  //}
  //return result;
  
  try {
    int size_out, size = meshIFace.getNumOfType( activeSet, TSTTM::EntityType_ALL_TYPES );
    sidl::array<void*> handles = alloc_sidl_vector<void*>(size);
    meshIFace.getEntities( activeSet, 
                           TSTTM::EntityType_ALL_TYPES, TSTTM::EntityTopology_ALL_TOPOLOGIES,
                           handles, size_out );
    if (size_out != size) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return 0;
    }
    return new SIDLIterator( handles );
  }
  catch (TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

//************ Vertex Properties ********************
// Returns true or false, indicating whether the vertex
// is allowed to be repositioned.  True indicates that the vertex
// is fixed and cannot be moved.  Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
bool MeshTSTTImpl::vertex_is_fixed( VertexHandle vertex, MsqError &err )
{
  try {
    return entIFace.getBoolData( vertex, fixedTag );
  }
  catch( TSTT::Error& tstt_err ) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return true;
  }
}

// Returns true or false, indicating whether the vertex
// is on the boundary.  Boundary nodes may be treated as
// a special case by some algorithms or culling methods.
// Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
void MeshTSTTImpl::vertices_are_on_boundary(
  VertexHandle vert_array[], bool bool_array[],
  size_t num_vtx, MsqError &err)
{
  try {
    int32_t lower=0, upper=num_vtx-1, stride=1;
    
    sidl::array<void*> vert_wrapper;
    vert_wrapper.borrow( vert_array, 1, &lower, &upper, &stride );
    
    /* SIDL's internal representation for an array<bool> is an 
       array of integers. Cannot use the input array directly
       (with array.borrow(..)), need to make a copy instead.

    sidl::array<bool> bool_wrapper;
    bool_wrapper.borrow( bool_array, 1, &lower, &upper, &stride );
    
    int num_bools = num_vtx;
    arrIFace.getBoolArrData( vert_wrapper, num_vtx, fixedTag, bool_array, num_bools );
    if (num_vtx != (unsigned)num_bools)
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
    */
    
    sidl::array<bool> bools = alloc_sidl_vector<bool>(num_vtx);
    for (size_t i = 0; i < num_vtx; ++i)
      bools.set( i, bool_array[i] );
    
    int num_bools = num_vtx;
    arrIFace.getBoolArrData( vert_wrapper, num_vtx, fixedTag, bools, num_bools );
    if (num_vtx != (unsigned)num_bools)
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

// Get vertex coordinates 
void MeshTSTTImpl::vertices_get_coordinates(
  Mesquite::Mesh::VertexHandle vert_array[],
  MsqVertex* const &coordinates, const size_t &num_vtx, MsqError &err)
{
  double* dbl_array = 0;
  
  try {
    int dim = meshIFace.getGeometricDim();
    int dbl_len = dim * num_vtx;

    int32_t lower=0, upper=num_vtx-1, stride=1;
    sidl::array<void*> vertex_wrapper;
    vertex_wrapper.borrow( vert_array, 1, &lower, &upper, &stride );
    dbl_array = new double[dbl_len];
    sidl::array<double> dbl_wrapper;
    upper = dbl_len = 1;
    dbl_wrapper.borrow( dbl_array, 1, &lower, &upper, &stride );
    
    TSTTM::StorageOrder order = meshIFace.getDfltStorage();
    meshIFace.getVtxArrCoords( vertex_wrapper, num_vtx, order, dbl_wrapper, dbl_len );
    if ((unsigned)dbl_len != dim*num_vtx) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      delete [] dbl_array;
      return;
    }
    
    if (dim == 2)
    {
      if (order == TSTTM::StorageOrder_INTERLEAVED)
      {
        double* iter = dbl_array;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].x(*iter); ++iter;
          coordinates[i].y(*iter); ++iter;
          coordinates[i].z(0);
        }
      }
      else
      {
        double *xiter = dbl_array;
        double *yiter = dbl_array + num_vtx;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].x(*xiter); ++xiter;
          coordinates[i].y(*yiter); ++yiter;
          coordinates[i].z(0);
        }
      }
    }
    else if(dim == 3)
    {
      if (order == TSTTM::StorageOrder_INTERLEAVED)
      {
        double* iter = dbl_array;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].set(iter);
          iter += 3;
        }
      }
      else
      {
        double *xiter = dbl_array;
        double *yiter = dbl_array + num_vtx;
        double *ziter = dbl_array + 2*num_vtx;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].x(*xiter); ++xiter;
          coordinates[i].y(*yiter); ++yiter;
          coordinates[i].z(*ziter); ++ziter;
        }
      }
    }
    else
    {
      MSQ_SETERR(err)(MsqError::INVALID_STATE, "TSTT database dimension = %d", dim);
      delete [] dbl_array;
      return;
    }
      
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
  
  delete [] dbl_array;
}

void MeshTSTTImpl::vertex_set_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  const Vector3D &coordinates, MsqError &err)
{
  try {
    int32_t lower = 0, upper = 2, stride = 1;
    sidl::array<double> coords;
    coords.borrow( const_cast<double*>(coordinates.to_array()), 
                   1, &lower, &upper, &stride );
    modIFace.setVtxCoords( vertex, coords );
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

// Each vertex has a byte-sized flag that can be used to store
// flags.  This byte's value is neither set nor used by the mesh
// implementation.  It is intended to be used by Mesquite algorithms.
// Until a vertex's byte has been explicitly set, its value is 0.
void MeshTSTTImpl::vertex_set_byte (
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char byte, MsqError &err)
{
  try {
    entIFace.setData( vertex, byteTag, &byte, 1 );
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::vertices_set_byte (
  VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
    // TSTT implementations seem to be inconsistant with
    // regard to setting opaque tags.  Set one at a time
    // to be safe.  This would be much easier if Mesquite
    // used a TSTT-defined type for the data, rather than
    // a single byte.
  try {
    VertexHandle *const end = vert_array + array_size;
    for ( ; vert_array != end; ++vert_array, ++byte_array)
      entIFace.setData( *vert_array, byteTag, byte_array, 1 );
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

// Retrieve the byte value for the specified vertex or vertices.
// The byte value is 0 if it has not yet been set via one of the
// *_set_byte() functions.
void MeshTSTTImpl::vertex_get_byte(
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char *byte, MsqError &err)
{
  try {
    int size;
    entIFace.getData( vertex, byteTag, (void*)byte, size );
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::vertices_get_byte(
  VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
    // TSTT implementations seem to be inconsistant with
    // regard to setting opaque tags.  Set one at a time
    // to be safe.  This would be much easier if Mesquite
    // used a TSTT-defined type for the data, rather than
    // a single byte.
  try {
    VertexHandle *const end = vert_array + array_size;
    int size;
    for ( ; vert_array != end; ++vert_array, ++byte_array)
      entIFace.getData( *vert_array, byteTag, (void*)byte_array, size );
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


//**************** Vertex Topology *****************


// Gets the number of elements attached to this vertex.
// Useful to determine how large the "elem_array" parameter
// of the vertex_get_attached_elements() function must be.
size_t MeshTSTTImpl::vertex_get_attached_element_count(
  VertexHandle vertex, MsqError &err)
{
  try {
    sidl::array<void*> junk;
    int face_size = 0, region_size = 0;
    entIFace.getEntAdj( vertex, TSTTM::EntityType_FACE, junk, face_size );
    entIFace.getEntAdj( vertex, TSTTM::EntityType_REGION, junk, region_size );
    return face_size + region_size;
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

// Gets the elements attached to this vertex.
void MeshTSTTImpl::vertex_get_attached_elements(
  VertexHandle vertex,
  ElementHandle* elem_array,
  size_t sizeof_elem_array, MsqError &err)
{
  try {
    sidl::array<void*> elem_wrapper;
    int face_size = 0, region_size = 0;
    entIFace.getEntAdj( vertex, TSTTM::EntityType_FACE, elem_wrapper, face_size );
    entIFace.getEntAdj( vertex, TSTTM::EntityType_REGION, elem_wrapper, region_size );
    
    if (sizeof_elem_array < (size_t)(face_size + region_size))
    {
      MSQ_SETERR(err)("Insufficient space in array", MsqError::OUT_OF_MEMORY);
      return;
    }
    
    int32_t lower = 0, upper, stride = 1;
    
    if (face_size)
    {
      upper = face_size - 1;
      elem_wrapper.borrow( elem_array, 1, &lower, &upper, &stride );
      entIFace.getEntAdj( vertex, TSTTM::EntityType_FACE, elem_wrapper, face_size );
    }
    
    if (region_size)
    {
      upper = region_size - 1;
      elem_wrapper.borrow( elem_array + face_size, 1, &lower, &upper, &stride );
      entIFace.getEntAdj( vertex, TSTTM::EntityType_REGION, elem_wrapper, region_size );
    }
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


// Gets the number of vertices in this element.
// This data can also be found by querying the
// element's topology and getting the number
// of vertices per element for that topology type.
size_t MeshTSTTImpl::element_get_attached_vertex_count(
  ElementHandle elem,
  MsqError &err)
{
  try {
    sidl::array<void*> junk;
    int result = 0;
    entIFace.getEntAdj( elem, TSTTM::EntityType_VERTEX, junk, result );
    return result;
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

/** Get connectivity
 *\param elements - Array of length num_elems containing elements
 *                  handles of elements for which connectivity is to
 *                  be queried.
 *\param vertices - Array of vertex handles in connectivity list.
 *\param offsets  - Indices into \ref indices array, one per element.
 *\param indices  - Indices into \ref vertex_handles
 */
void MeshTSTTImpl::elements_get_attached_vertices(
  ElementHandle *elements,
  size_t num_elems,
  VertexHandle *vertices,
  size_t &vert_len,
  size_t *indices,
  size_t &indices_len,
  size_t *offsets,
  Mesquite::MsqError &err)
{
  if (num_elems == 0)
    return;

  try {
    int32_t lower = 0, upper, stride = 1;
    
    int vert_count = vert_len;
    int index_count = indices_len;
    int off_count = num_elems;
    sidl::array<void*> vert_wrapper;
    sidl::array<int> index_wrapper, offset_wrapper;
    upper = vert_len - 1;
    vert_wrapper.borrow( vertices, 1, &lower, &upper, &stride );
    
    // If sizeof(size_t) == sizeof(int), can avoid copying
    if (sizeof(size_t) == sizeof(int))
    {
      upper = indices_len - 1;
      index_wrapper.borrow( (int*)indices, 1, &lower, &upper, &stride );
      upper = num_elems - 1;
      offset_wrapper.borrow( (int*)offsets, 1, &lower, &upper, &stride );
    }
    else
    {
      index_wrapper = alloc_sidl_vector<int>(indices_len);
      offset_wrapper = alloc_sidl_vector<int>(num_elems);
    } 

    // The only TSTT API that I can find that gives the data in the
    // format Mesquite wants takes an entity set as input.  Create
    // entity set from handle array.
    void* set;
    setIFace.createEntSet( true, set );
    sidl::array<void*> elem_wrapper;
    upper = num_elems-1;
    elem_wrapper.borrow( elements, 1, &lower, &upper, &stride );
    setIFace.addEntArrToSet( set, elem_wrapper, num_elems );
    
    // Get the data.  
    meshIFace.getAdjEntities( set, 
                              TSTTM::EntityType_ALL_TYPES, 
                              TSTTM::EntityTopology_ALL_TOPOLOGIES,
                              TSTTM::EntityType_VERTEX,
                              vert_wrapper, vert_count,
                              offset_wrapper, off_count,
                              index_wrapper, index_count );
    
    // Destroy the temporary entity set
    setIFace.destroyEntSet( set );
    
    // Check sizes that were passed back from TSTT database
    if ((size_t)off_count != num_elems)
    {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return;
    }
    if ((size_t)vert_count > vert_len || (size_t)index_count > indices_len)
    {
      // Pass required sizes back
      vert_len = vert_count;
      indices_len = index_count;
      MSQ_SETERR(err)("Insufficient space in array", MsqError::OUT_OF_MEMORY) ;
      return;
    }
    
    // Pass back the actual length of the data
    vert_len = vert_count;
    indices_len = index_count;
    
    // if sizes not same, need to cast/copy data.
    if (sizeof(size_t) != sizeof(int))
    {
      copy_from_sidl( index_wrapper, indices );
      copy_from_sidl( offset_wrapper, offsets );
    }
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::get_all_sizes( size_t& vertex_count,
                                  size_t& element_count,
                                  size_t& vertex_use_count,
                                  MsqError& err )
{
  try {
    
    sidl::array<void*> handles;
    sidl::array<int> offsets;
    sidl::array<int> indices;
    int num_handles, num_offsets, num_indices;
    
    meshIFace.getAdjEntities( activeSet, 
                              TSTTM::EntityType_ALL_TYPES,
                              TSTTM::EntityTopology_ALL_TOPOLOGIES,
                              TSTTM::EntityType_VERTEX,
                              handles, num_handles,
                              offsets, num_offsets,
                              indices, num_indices );
    
    vertex_count = num_handles;
    element_count = num_offsets-1;
    vertex_use_count = num_indices;
    
    if (element_count != (size_t)meshIFace.getNumOfType(activeSet, TSTTM::EntityType_ALL_TYPES))
    {
      MSQ_SETERR(err)( MsqError::INTERNAL_ERROR );
    }
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::get_all_mesh( VertexHandle* vert_array,  size_t vert_len,
                                 ElementHandle* elem_array, size_t elem_len,
                                 size_t* offset_array,      size_t offset_len,
                                 size_t* conn_array,        size_t conn_len ,
                                 MsqError& err )
{
  try {
    int32_t lower = 0, upper, stride = 1;
    sidl::array<void*> vertices;
    sidl::array<void*> elements;
    sidl::array<int>   offsets;
    sidl::array<int>   indices;
    
    upper = vert_len;
    vertices.borrow( vert_array, 1, &lower, &upper, &stride );
    upper = elem_len;
    elements.borrow( elem_array, 1, &lower, &upper, &stride );
    
    int num_elem_out;
    meshIFace.getEntities( activeSet, 
                           TSTTM::EntityType_ALL_TYPES,
                           TSTTM::EntityTopology_ALL_TOPOLOGIES,
                           elements,
                           num_elem_out );
    if ((unsigned)num_elem_out > elem_len)
    {
      MSQ_SETERR(err)("Insufficient space in array", MsqError::OUT_OF_MEMORY);
      return;
    }
    
    if (sizeof(int) == sizeof(size_t))
    {
      upper = offset_len;
      offsets.borrow( (int*)offset_array, 1, &lower, &upper, &stride );
      upper = conn_len;
      indices.borrow( (int*)conn_array, 1, &lower, &upper, &stride );
    }
    else
    {
      offsets = alloc_sidl_vector<int>(offset_len);
      indices = alloc_sidl_vector<int>(conn_len);
    }
    
    int num_vert_out, num_uses_out;
    meshIFace.getAdjEntities( activeSet,
                              TSTTM::EntityType_ALL_TYPES,
                              TSTTM::EntityTopology_ALL_TOPOLOGIES,
                              TSTTM::EntityType_VERTEX,
                              vertices, num_vert_out,
                              offsets,  num_elem_out,
                              indices,  num_uses_out );
    
    if ((unsigned)num_vert_out > vert_len ||
        (unsigned)num_elem_out > offset_len ||
        (unsigned)num_uses_out > conn_len) 
    {
      MSQ_SETERR(err)("Insufficient space in array", MsqError::OUT_OF_MEMORY);
      return;
    }
    
    if (sizeof(int) != sizeof(size_t))
    {
      copy_from_sidl( offsets, offset_array );
      copy_from_sidl( indices, conn_array );
    }
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}
      

// Returns the topology of the given entity.
EntityTopology MeshTSTTImpl::element_get_topology(
  ElementHandle element, MsqError &err)
{
  try {
    TSTTM::EntityTopology topo = entIFace.getEntTopo( element );
    return topologyMap[topo];
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return MIXED;
  }
}

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
void MeshTSTTImpl::elements_get_topologies(
  ElementHandle *element_handle_array,
  EntityTopology *element_topologies,
  size_t num_elements, MsqError &err)
{
  try {
    assert(sizeof(EntityTopology) == sizeof(TSTTM_EntityTopology__enum));
    
    int32_t lower = 0, upper = num_elements, stride = 1;
    sidl::array<TSTTM::EntityTopology> topologies;
    sidl::array<void*> handles;
    handles.borrow( element_handle_array, 1, &lower, &upper, &stride );
    topologies.borrow( reinterpret_cast<TSTTM_EntityTopology__enum*>(element_topologies), 
                       1, &lower, &upper, &stride );
    int topo_len_out;
    arrIFace.getEntArrTopo( handles, num_elements, topologies, topo_len_out );
    if ((size_t)topo_len_out != num_elements) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return;
    }
    
    EntityTopology *iter = element_topologies, 
              *const end = element_topologies + num_elements;
    for (; iter != end; ++iter)
      *iter = topologyMap[*iter];
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void MeshTSTTImpl::release_entity_handles(
  Mesquite::Mesh::EntityHandle */*handle_array*/,
  size_t /*num_handles*/, MsqError &/*err*/)
{
    // Do nothing
}

// Instead of deleting a Mesh when you think you are done,
// call release().  In simple cases, the implementation could
// just call the destructor.  More sophisticated implementations
// may want to keep the Mesh object to live longer than Mesquite
// is using it.
void MeshTSTTImpl::release()
{
}

//**************** Tags ****************
TagHandle MeshTSTTImpl::tag_create( const msq_std::string& name, 
                                    TagType type, unsigned length,
                                    const void* default_val,
                                    MsqError& err )
{
  TSTT::TagValueType tstt_type;
  size_t size = 0;
  switch (type) {
    case BYTE:   size = sizeof(char  ); tstt_type = TSTT::TagValueType_OPAQUE;        break;
    case BOOL:   size = sizeof(bool  ); tstt_type = TSTT::TagValueType_BOOLEAN;       break;
    case INT:    size = sizeof(int   ); tstt_type = TSTT::TagValueType_INTEGER;       break;
    case DOUBLE: size = sizeof(double); tstt_type = TSTT::TagValueType_DOUBLE;        break;
    case HANDLE: size = sizeof(void* ); tstt_type = TSTT::TagValueType_ENTITY_HANDLE; break;
    default:
      MSQ_SETERR(err)("Invalid tag type", MsqError::INVALID_ARG );
      return 0;
  }
  size *= length;
  
  try {
    void* handle = 0;
    entIFace.createTag( name, size, tstt_type, const_cast<void*>(default_val), handle );
    return handle;
  } 
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

void MeshTSTTImpl::tag_delete( TagHandle handle, MsqError& err )
{
  try {
    entIFace.destroyTag( handle, true );
  } 
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

TagHandle MeshTSTTImpl::tag_get( const msq_std::string& name, MsqError& err )
{
  try {
    return entIFace.getTagHandle( name );
  } 
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

void MeshTSTTImpl::tag_properties( TagHandle handle,
                                   msq_std::string& name,
                                   TagType& type_out,
                                   unsigned& length_out,
                                   MsqError& err )
{
  size_t size;
  TSTT::TagValueType type;
  try {
    name = entIFace.getTagName( handle );
    size = entIFace.getTagSize( handle );
    type = entIFace.getTagType( handle );
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
  
  size_t tsize;
  switch (type) {
    case TSTT::TagValueType_OPAQUE       : tsize = sizeof(char  ); type_out = BYTE  ; break;
    case TSTT::TagValueType_BOOLEAN      : tsize = sizeof(bool  ); type_out = BOOL  ; break;
    case TSTT::TagValueType_INTEGER      : tsize = sizeof(int   ); type_out = INT   ; break;
    case TSTT::TagValueType_DOUBLE       : tsize = sizeof(double); type_out = DOUBLE; break;
    case TSTT::TagValueType_ENTITY_HANDLE: tsize = sizeof(void* ); type_out = HANDLE; break;
    default:
      MSQ_SETERR(err)("Unsupported TSTT tag type", MsqError::NOT_IMPLEMENTED );
      return;
  }
  
  if (size % tsize != 0)
  {
    MSQ_SETERR(err)("Invalid TSTT tag size", MsqError::INTERNAL_ERROR );
    return;
  }
  
  length_out = size / tsize;
}

void MeshTSTTImpl::tag_set_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         const void* data,
                                         MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}

void MeshTSTTImpl::tag_set_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        const void* data,
                                        MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}
    
void MeshTSTTImpl::tag_set_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 const void* data,
                                 MsqError& err )
{
  try {
    size_t size = entIFace.getTagSize( tag );
    sidl::array<void*> values = alloc_sidl_vector<void*>(num_elems);
    const char* ptr = (const char*)data;
    for (size_t i = 0; i < num_elems; ++i, ptr += size)
      values.set( i, (void*)ptr );
    
    int32_t lower = 0, upper = num_elems - 1, stride = 1;
    sidl::array<void*> handles;
    handles.borrow( const_cast<void**>(array), 1, &lower, &upper, &stride );
    
    arrIFace.setArrData( handles, num_elems, tag, values, num_elems, size );
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


void MeshTSTTImpl::tag_get_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         void* data,
                                         MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}

void MeshTSTTImpl::tag_get_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        void* data,
                                        MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}
    
void MeshTSTTImpl::tag_get_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 void* data,
                                 MsqError& err )
{
  try {
    size_t size = entIFace.getTagSize( tag );
    sidl::array<void*> values = alloc_sidl_vector<void*>(num_elems);
    char* ptr = reinterpret_cast<char*>(data);
    for (size_t i = 0; i < num_elems; ++i, ptr += size)
      values.set( i, (void*)ptr );
    
    int32_t lower = 0, upper = num_elems - 1, stride = 1;
    sidl::array<void*> handles;
    handles.borrow( const_cast<void**>(array), 1, &lower, &upper, &stride );
    
    int num_vals, size_out;
    arrIFace.getArrData( handles, num_elems, tag, values, num_vals, size_out );
    if ((size_t)num_vals != num_elems || size != (size_t)size_out)
    {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
    }
  }
  catch(::TSTT::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

} // namespace Mesquite
