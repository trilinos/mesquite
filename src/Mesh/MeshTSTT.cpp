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

#include <sidl_cxx.hh>
#include "TSTTB.hh"
#include "TSTTM.hh"
#include "MeshTSTT.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "MeshInterface.hpp"
#include "TSTTUtil.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <algorithm.h>
#else
# include <algorithm>
#endif


#define OPAQUE_TYPE_OPAQUE_PADDED 0
#define OPAQUE_TYPE_OPAQUE_PACKED 1
#define OPAQUE_TYPE_CHAR          2
#define OPAQUE_TYPE_UCHAR         3
#define OPAQUE_TYPE_BYTE          OPAQUE_TYPE_UCHAR
#define TSTT_OPAQUE_TAG_TYPE OPAQUE_TYPE_CHAR


namespace Mesquite
{



/*************************************************************************
 *                          Iterator Definitions
 *
 * This file contains three alternate iterator implementations.  After
 * some performance testing, two should be removed.
 ************************************************************************/

/** \brief Wrapper around single-entity TSTT interator */
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
  catch (TSTTB::Error& tstt_err) {
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
  catch (TSTTB::Error& tstt_err) {
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
  catch (TSTTB::Error& tstt_err) {
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

/** \brief TSTT iterator using array-iterator interface and buffer of handles */
class TSTTArrIter : public EntityIterator
{
  private:
  
    sidl::array<void*> handleArray;
    void* tsttIter;
    int index, count;
    bool notAtEnd;
    TSTTM::Arr& myMesh;
    
    inline void get_next_array() 
    {
      index = count = 0;
      notAtEnd = myMesh.getNextEntArrIter( tsttIter, handleArray, count ) && count;
    }
    
  public:
  
    TSTTArrIter( TSTTM::Arr& mesh, 
                 void* meshset, 
                 TSTTM::EntityType type,
                 TSTTM::EntityTopology topo,
                 unsigned buffer_count = 1024 );
    
      /**\brief reset iterator */
    virtual void restart();
      /**\brief get current entity handle */
    virtual Mesh::EntityHandle operator*() const;
      /**\brief check if any remaining entity handles */
    virtual bool is_at_end() const;
      /**\biref step */
    virtual void operator++();
    
    virtual ~TSTTArrIter();
};

TSTTArrIter::TSTTArrIter( TSTTM::Arr& mesh, 
                          void* meshset, 
                          TSTTM::EntityType type,
                          TSTTM::EntityTopology topo,
                          unsigned buffer_count )
                        : handleArray( alloc_sidl_vector<void*>(buffer_count) ),
                          tsttIter(0),
                          index(0),
                          count(0),
                          notAtEnd(false),
                          myMesh( mesh )
{
  mesh.initEntArrIter( meshset, type, topo, buffer_count, tsttIter );
  get_next_array();
}

TSTTArrIter::~TSTTArrIter() 
{
  try {
    if (tsttIter)
      myMesh.endEntArrIter( tsttIter );
  }
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}

void TSTTArrIter::restart()
{
  try {
    myMesh.resetEntArrIter( tsttIter );
    get_next_array();
  }
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}

Mesh::EntityHandle TSTTArrIter::operator*() const
{
  return handleArray.get(index);
}

bool TSTTArrIter::is_at_end() const
{
  return index == count && !notAtEnd;
}

void TSTTArrIter::operator++()
{
  try {
    ++index;
    if (index == count && notAtEnd)
      get_next_array();
  }                        
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}
     
      



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
       *                  mesh elements and vertices for which quality 
       *                  is to be improved.
       */
    virtual void set_active_set( void* element_set, MsqError& );
    

      /**\brief Get dimension of vertex coordinates (2D vs. 3D). */
    virtual int get_geometric_dimension(Mesquite::MsqError &/*err*/);
    
    /** \brief Get handles for all elemnents */
    virtual void get_all_elements( msq_std::vector<ElementHandle>& elements, 
                                   MsqError& err );
    
    /** \brief Get handles for all vertices */
    virtual void get_all_vertices( msq_std::vector<VertexHandle>& vertices, 
                                   MsqError& err );
    
      /**\brief Create iterator for vertices in active set */
    virtual VertexIterator* vertex_iterator(MsqError &err);
    
      /**\brief Create iterator for elements in active set */
    virtual ElementIterator* element_iterator(MsqError &err);

      /**\brief Query "fixed" flag for a vertex */
    virtual void vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                          bool fixed_flag_array[],
                                          size_t num_vtx, 
                                          MsqError &err);
 
      /**\brief Get vertex coordinates */
    virtual void vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx, 
                                           MsqError &err);
      /**\brief Set vertex coordinates */
    virtual void vertex_set_coordinates( VertexHandle vertex,
                                         const Vector3D &coordinates, 
                                         MsqError &err);
    
      /**\brief Set vertex mark */
    virtual void vertex_set_byte( VertexHandle vertex,
                                  unsigned char byte, 
                                  MsqError &err);
      /**\brief Set vertex mark */
    virtual void vertices_set_byte( const VertexHandle *vert_array,
                                    const unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err);
    
      /**\brief Get vertex mark */
    virtual void vertex_get_byte( VertexHandle vertex,
                                  unsigned char *byte, 
                                  MsqError &err);
      /**\brief Get vertex mark */
    virtual void vertices_get_byte( const VertexHandle *vert_array,
                                    unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err);
    
      /**\brief Get vertex-to-element adjacencies */
    virtual void vertices_get_attached_elements( const VertexHandle* vertex_array,
                                                 size_t num_vertices,
                                                 msq_std::vector<ElementHandle>& elements,
                                                 msq_std::vector<size_t>& offsets,
                                                 MsqError& err );
    
      /**\brief Get element connectivity */
    virtual void elements_get_attached_vertices( const ElementHandle *elem_handles,
                                                 size_t num_elems,
                                                 msq_std::vector<VertexHandle>& vertices,
                                                 msq_std::vector<size_t>& offsets,
                                                 MsqError& err );
    
  
      /**\brief Return topology type enum for an array of elements */
    virtual void elements_get_topologies( const ElementHandle *element_handle_array,
                                          EntityTopology *element_topologies,
                                          size_t num_elements, 
                                          MsqError &err );
    
//**************** Memory Management ****************
      /**\brief no-op */ 
    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles, 
                                         MsqError &err );
    
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
                                       
  protected:
        
    void set_int_tag( void* tag, void* meshset, int value, MsqError& err );

    /** Populate \ref inputElements from \ref elementSet */
    void populate_input_elements( MsqError& err );

      /** \brief  Call TSTTM::Arr::getEntArrAdj
       *
       * Common code for \ref vertices_get_attached_elements and 
       * \ref elements_get_attached_vertices
       *
       *\param source      Array of handles of source entities to query from
       *\param num_source  The length of \ref source
       *\param target_type The type of entity to query for
       *\param target      The output list of adjacent entities
       *\param offsets     For each entity in \ref source, the offset in 
       *                   \ref target at which the corresponding adjacent
       *                   entities are stored. (output)
       */
    void get_adjacent_entities( void* const* source,
                                size_t num_source,
                                TSTTM::EntityType target_type,
                                msq_std::vector<void*>& target,
                                msq_std::vector<size_t>& offsets,
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
    TSTTM::Entity entIFace;
    TSTTB::EntTag tagIFace;
    /** TSTT interface for multi-entity (array) queries */
    TSTTM::Arr arrIFace;
    TSTTB::ArrTag arrTagIFace;
    /** TSTT interface for modifying mesh */
    TSTTM::Modify modIFace;
    /** TSTT interface for entity set operations */
    TSTTB::EntSet setIFace;
    
    /** Have mesh */
    bool haveMesh;
    /** TSTTM entity set handle for elements to improve */
    void* elementSet;
    /** TSTTM entity set handle for nodes to move */
    void* nodeSet;
    /** std::set containing elements in elementSet, used
     *  to constrain vertex->element adjaceny queries to
     *  only those elements that are in the input element set.
     */
    msq_std::vector<void*> inputElements;
    
    /** The type of elements contained in the input element set.
     * Should be one of:
     * - TSTTM::EntityType_REGION    - volume elements
     * - TSTTM::EntityType_FACE      - face/2d elements
     * - TSTTM::EntityType_ALL_TYPES - mixed volume and face elements
     */
    TSTTM::EntityType inputSetType;
    
    /** Handle for tag used to hold vertex byte */
    TagHandle byteTag; 
    /** Tag was created in constructor */
    bool createdByteTag;
    /** Handle for tag used to hold vertex-fixed flag */
    TagHandle fixedTag;
    /** Fixed tag was created in constructor */
    bool createdFixedTag;
    /** Handle for the tag used internally to remove duplicates from lists */
//    TagHandle vertexIndexTag;
    /** vertexIndexTag was created in constructor */
//    bool createdVertexIndexTag;
    
    /** Map TSTTM::EntityTopology to Mesquite::EntityTopology */
    EntityTopology topologyMap[TSTTM::EntityTopology_ALL_TOPOLOGIES+1];
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
    return 0;
  }
  result->set_active_set( meshset, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  return result;
}

MeshTSTT* MeshTSTT::create( TSTTM::Mesh& mesh, MsqError& err )
{
  MeshTSTT* result = new MeshTSTTImpl( mesh, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  return result;
}

MeshTSTT::~MeshTSTT() {}


MeshTSTTImpl::MeshTSTTImpl(TSTTM::Mesh& tstt_mesh, Mesquite::MsqError& err) 
  : meshIFace(tstt_mesh), 
    elementSet(0), nodeSet(0), 
    inputSetType( TSTTM::EntityType_ALL_TYPES ),
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
      // Get tag interface
    tagIFace = tstt_mesh;
    if (!tagIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTB::EntTag",
                       MsqError::INVALID_STATE );
      return;
    }
      // Get tag interface
    arrTagIFace = tstt_mesh;
    if (!arrTagIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTB::ArrTag",
                       MsqError::INVALID_STATE );
      return;
    }
    
      // Get tag for fixed flag
    try {
      fixedTag = tagIFace.getTagHandle( VERTEX_FIXED_TAG_NAME );

        // Double-check types incase tag already existed
      if (tagIFace.getTagSizeBytes(fixedTag) != sizeof(int) || 
          tagIFace.getTagType(fixedTag) != TSTTB::TagValueType_INTEGER)
      {
        MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                         "Tag \"%s\" exists with invalid type/size", 
                         VERTEX_FIXED_TAG_NAME );
        return;
      }
    } catch (...) { 
      fixedTag = 0;
    }
    
      // Get/create tag for vertex byte
    try {
      byteTag = tagIFace.getTagHandle( VERTEX_BYTE_TAG_NAME );
    } catch (...) {}

    if (!byteTag) {
      tagIFace.createTag( VERTEX_BYTE_TAG_NAME, 1, TSTTB::TagValueType_INTEGER, byteTag );
      createdByteTag = true;
    } 
      // Double-check types incase tag already existed
    if (tagIFace.getTagSizeBytes(byteTag) != sizeof(int) || 
        tagIFace.getTagType(byteTag) != TSTTB::TagValueType_INTEGER)
    {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                       "Tag \"%s\" exists with invalid type/size", 
                       VERTEX_BYTE_TAG_NAME );
      return;
    }
    
      // Clear vertex byte tag
    set_int_tag( byteTag, nodeSet, 0, err );
    if (err)
      return;
    
      // Get/create tag for vertex index
/*
    //try {
    //  vertexIndexTag = tagIFace.getTagHandle( VERTEX_INDEX_TAG_NAME );
    //} catch (...) {}

    if (!vertexIndexTag) {
      tagIFace.createTag( VERTEX_INDEX_TAG_NAME, sizeof(int), TSTTB::TagValueType_INTEGER, byteTag );
      createdVertexIndexTag = true;
    } 
      // Double-check types incase tag already existed
    if (tagIFace.getTagSize(vertexIndexTag) != sizeof(int) || 
        tagIFace.getTagType(vertexIndexTag) != TSTTB::TagValueType_INTEGER)
    {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                       "Tag \"%s\" exists with invalid type/size", 
                       VERTEX_INDEX_TAG_NAME );
      return;
    }
*/
  }
  catch (TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    try {
      if (createdFixedTag)
        tagIFace.destroyTag( fixedTag, false );
      if (createdByteTag)
        tagIFace.destroyTag( byteTag, false );
//      if (createdVertexIndexTag)
//        tagIFace.destroyTag( vertexIndexTag, false );
    } 
    catch (...) {}
  }
  catch(...) {
    MSQ_SETERR(err)("Uknown exception",MsqError::INTERNAL_ERROR);
  }
}

Mesquite::MeshTSTTImpl::~MeshTSTTImpl() 
{
  try {
    if (elementSet) 
      setIFace.destroyEntSet( elementSet );
    if (nodeSet)
      setIFace.destroyEntSet( nodeSet );
      
    if (createdFixedTag)
      tagIFace.destroyTag( fixedTag, false );
    if (createdByteTag)
      tagIFace.destroyTag( byteTag, false );
//    if (createdVertexIndexTag)
//      tagIFace.destroyTag( vertexIndexTag, false );
  } 
  catch (TSTTB::Error& tstt_err ) {
    process_tstt_error(tstt_err);
    throw;
  }
}


void MeshTSTTImpl::set_int_tag( void* tag,
                                void* elem_set, 
                                int value, 
                                MsqError& err )
{
  const unsigned BUFFER_COUNT = 1024;
  
  sidl::array<int> value_array( alloc_sidl_vector<int>( BUFFER_COUNT, value ) );
  
  sidl::array<void*> handle_array( alloc_sidl_vector<void*>( BUFFER_COUNT ) );
  int count = BUFFER_COUNT;
  void* iter = 0;
  bool more;
  
  try {
    
    arrIFace.initEntArrIter( elem_set, 
                             TSTTM::EntityType_VERTEX, 
                             TSTTM::EntityTopology_POINT, 
                             BUFFER_COUNT, iter );
                             
    do {
      more = arrIFace.getNextEntArrIter( iter, handle_array, count );
      if (count > 0)
        arrTagIFace.setIntArrData( handle_array, count, tag, value_array, count );
    } while (more);
    
    arrIFace.endEntArrIter( iter );
  }
  
  catch (TSTTB::Error &tstt_err) {
    if( iter ) try {  arrIFace.endEntArrIter( iter ); } catch (...) {}
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


void MeshTSTTImpl::set_active_set( void* elem_set, MsqError& err )
{
  const int ELEM_BUFFER_SIZE = 1024;
  const int NODE_BUFFER_SIZE = 27 * ELEM_BUFFER_SIZE; 
  sidl::array<void*> elements( alloc_sidl_vector<void*>( ELEM_BUFFER_SIZE ) );
  sidl::array<void*>    nodes( alloc_sidl_vector<void*>( NODE_BUFFER_SIZE ) );
  sidl::array<int>    offsets( alloc_sidl_vector<int  >( ELEM_BUFFER_SIZE+1 ) );
  void* iter = 0;
 
  try {
    
      // Create new sets for nodes and elements
    if (elementSet)
    {
      setIFace.destroyEntSet( elementSet );
      elementSet = 0;
    }
    if (nodeSet)
    {
      setIFace.destroyEntSet( nodeSet );
      nodeSet = 0;
    }
    setIFace.createEntSet( false, elementSet );
    setIFace.createEntSet( false, nodeSet );
    
      // Iterate over set twice, once for FACEs and once for REGIONs
    bool have_faces = false, have_regions = false;
    for (int i = 0; i < 2; ++i)
    {
      TSTTM::EntityType type = i ? TSTTM::EntityType_REGION : 
                                   TSTTM::EntityType_FACE;
      bool& have_some = i ? have_regions : have_faces;
                                   
      arrIFace.initEntArrIter( elem_set, 
                               type, 
                               TSTTM::EntityTopology_ALL_TOPOLOGIES,
                               ELEM_BUFFER_SIZE, 
                               iter );
      
      int count = 0;
      bool more = false;
      do {
          // Add elements to element set
        more = arrIFace.getNextEntArrIter( iter, elements, count );
        if (!count) break;
        setIFace.addEntArrToSet( elements, count, elementSet );
        
          // Add nodes to node set
        int num_nodes, num_offsets;
        arrIFace.getEntArrAdj( elements, count, TSTTM::EntityType_VERTEX,
                               nodes, num_nodes, offsets, num_offsets );
        setIFace.addEntArrToSet( nodes, num_nodes, nodeSet );
        
        have_some = true;
      } while (more);
        
      arrIFace.endEntArrIter( iter );
      iter = 0;
      
    } // for (type)
    
    if (!have_faces)
      inputSetType = TSTTM::EntityType_REGION;
    else if (!have_regions)
      inputSetType = TSTTM::EntityType_FACE;
    else
      inputSetType = TSTTM::EntityType_ALL_TYPES;
    
    //set_fixed_tag( nodeSet, false, err ); MSQ_ERRRTN(err);
  }
  catch (TSTTB::Error &tstt_err) {
      // If an error occured, try to clean up anything created above

    try {
      if (iter) {
        arrIFace.endEntArrIter( iter );
        iter = 0;
      }
    }
    catch( ... ) {}
    
    try {
      if (elementSet) {
        setIFace.destroyEntSet( elementSet );
        elementSet = 0;
      }
    }
    catch( ... ) {}
    
    try {
      if (nodeSet) {
        setIFace.destroyEntSet( nodeSet );
        nodeSet = 0;
      }
    }
    catch( ... ) {}
    
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
  
    // clear cached data
  inputElements.clear();
}

void MeshTSTTImpl::populate_input_elements( MsqError& err ) 
{
  inputElements.clear();
  get_all_elements( inputElements, err );
  if (MSQ_CHKERR(err))
  {
    inputElements.clear();
    return;
  }
  
  msq_std::sort( inputElements.begin(), inputElements.end() );
} 

  

// Returns whether this mesh lies in a 2D or 3D coordinate system.
int MeshTSTTImpl::get_geometric_dimension(Mesquite::MsqError &err)
{
  try {
    return meshIFace.getGeometricDim();
  }
  catch (TSTTB::Error &tstt_err) {
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
    return new TSTTArrIter( arrIFace, 
                            nodeSet, 
                            TSTTM::EntityType_ALL_TYPES,
                            TSTTM::EntityTopology_ALL_TOPOLOGIES );
  }
  catch (TSTTB::Error &tstt_err) {
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
  try {
    return new TSTTArrIter( arrIFace, 
                            elementSet, 
                            TSTTM::EntityType_ALL_TYPES, 
                            TSTTM::EntityTopology_ALL_TOPOLOGIES );
  }
  catch (TSTTB::Error &tstt_err) {
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
void MeshTSTTImpl::vertices_get_fixed_flag(
  const VertexHandle vert_array[], 
  bool bool_array[],
  size_t num_vtx, MsqError &err)
{
    // If mesh does not contain a fixed tag, assume no vertices are fixed
  if (!fixedTag) {
    memset( bool_array, 0, num_vtx * sizeof(bool) );
    return;
  }
  
    // Get per-vertex flags from fixedTag
  try {
    
    if (num_vtx == 1)
    {
      bool_array[0] = (bool)tagIFace.getIntData( vert_array[0], fixedTag );
    }
    else
    {
      sidl::array<void*> vert_wrapper( convert_to_sidl_vector( const_cast<void**>(vert_array), num_vtx ) );
      sidl::array<int> bools = alloc_sidl_vector<int>(num_vtx);
    
      int num_bools = num_vtx;
      arrTagIFace.getIntArrData( vert_wrapper, num_vtx, fixedTag, bools, num_bools );
      if (num_vtx != (unsigned)num_bools)
        MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);

      for (size_t i = 0; i < num_vtx; ++i)
        bool_array[i] = bools.get(i);
    }
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

// Get vertex coordinates 
void MeshTSTTImpl::vertices_get_coordinates(
  const Mesquite::Mesh::VertexHandle vert_array[],
  MsqVertex* coordinates, 
  size_t num_vtx, 
  MsqError &err)
{
  double* dbl_array = 0;
  
  try {
    int dim = meshIFace.getGeometricDim();
    int dbl_len = dim * num_vtx;

    sidl::array<void*> vertex_wrapper( 
      convert_to_sidl_vector( const_cast<void**>(vert_array), num_vtx ) );
    dbl_array = new double[dbl_len];
    sidl::array<double> dbl_wrapper( convert_to_sidl_vector( dbl_array, dbl_len ));
    
    TSTTM::StorageOrder order = TSTTM::StorageOrder_UNDETERMINED;
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
      MSQ_SETERR(err)(MsqError::INVALID_STATE, "TSTTB database dimension = %d", dim);
      delete [] dbl_array;
      return;
    }
      
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
  
  delete [] dbl_array;
}

void MeshTSTTImpl::vertex_set_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  const Vector3D &coords, MsqError &err)
{
  try {
    modIFace.setVtxCoords( vertex, coords[0], coords[1], coords[2] );
  }
  catch(::TSTTB::Error &tstt_err) {
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
    int value = byte;
    tagIFace.setIntData( vertex, byteTag, value );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::vertices_set_byte (
  const VertexHandle *vert_array,
  const unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
    // TSTT implementations seem to be inconsistant with
    // regard to setting opaque tags.  Set one at a time
    // to be safe.  This would be much easier if Mesquite
    // used a TSTT-defined type for the data, rather than
    // a single byte.
  try {
    sidl::array<void*> handles( convert_to_sidl_vector( const_cast<void**>(vert_array), array_size ) );
    sidl::array<int> data(alloc_sidl_vector<int>(array_size));
    for (size_t i = 0; i < array_size; ++i)
      data.set( i, byte_array[i] );
    arrTagIFace.setIntArrData( handles, array_size, byteTag, data, array_size );
  }
  catch(::TSTTB::Error &tstt_err) {
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
    *byte = (unsigned char)tagIFace.getIntData( vertex, byteTag );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::vertices_get_byte(
  const VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
    // TSTT implementations seem to be inconsistant with
    // regard to setting opaque tags.  Set one at a time
    // to be safe.  This would be much easier if Mesquite
    // used a TSTT-defined type for the data, rather than
    // a single byte.
  try {
    sidl::array<void*> handles( convert_to_sidl_vector( const_cast<void**>(vert_array), array_size ));
    sidl::array<int> data( alloc_sidl_vector<int>(array_size) );
    int32_t junk;
    arrTagIFace.getIntArrData( handles, array_size, byteTag, data, junk );

    for (size_t i = 0; i< array_size; ++i )
      byte_array[i] = (unsigned char)data.get(i);
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


//**************** Topology *****************

void MeshTSTTImpl::get_adjacent_entities( void* const* source,
                                          size_t num_source,
                                          TSTTM::EntityType target_type,
                                          msq_std::vector<void*>& target,
                                          msq_std::vector<size_t>& offsets,
                                          MsqError& err )
{
  if (num_source == 0)
    return;
  
  int num_adj = 0, num_offset;
  sidl::array<int> sidl_offsets;
  
  sidl::array<void*> sidl_source = convert_to_sidl_vector( const_cast<void**>(source), num_source );
  offsets.resize( num_source + 1 );
  if (sizeof(size_t) == sizeof(int)) // avoid copy if possible
    sidl_offsets = convert_to_sidl_vector( reinterpret_cast<int*>(&offsets[0]), offsets.size() );
  
  bool have_adj = false;
    // If passed vector has allocated storage, try to use existing space
  if (target.capacity() >= num_source)
  {
    target.resize( target.capacity() );
    sidl::array<void*> sidl_target = convert_to_sidl_vector( &target[0], target.size() );

    try {
      arrIFace.getEntArrAdj( sidl_source, num_source,
                             target_type,
                             sidl_target, num_adj,
                             sidl_offsets, num_offset );
      assert( (unsigned)num_adj <= target.size() );
      target.resize( num_adj );
      have_adj = true;
    } catch (...) { }
  }
  
    // If implementation passed back a size, try that
  if (!have_adj && num_adj && (unsigned)num_adj > target.capacity())
  {
    target.resize( target.capacity() );
    sidl::array<void*> sidl_target = convert_to_sidl_vector( &target[0], target.size() );

    try {
      arrIFace.getEntArrAdj( sidl_source, num_source,
                             target_type,
                             sidl_target, num_adj,
                             sidl_offsets, num_offset );
      assert( (unsigned)num_adj <= target.size() );
      target.resize( num_adj );
      have_adj = true;
    } catch (...) { }
  }

    // Try with empty sidl array, and copy into elements vector
  if (!have_adj)
  {
    try {
      sidl::array<void*> sidl_target;
      arrIFace.getEntArrAdj( sidl_source, num_source,
                             target_type,
                             sidl_target, num_adj,
                             sidl_offsets, num_offset );

      target.resize( num_adj );
      copy_from_sidl( sidl_target, &target[0] );
    }
    catch(::TSTTB::Error &tstt_err) {
      MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
      return;
    }
  }
  
  // TSTT implementations seem to be inconsistent with regard to 
  // placing the last value on this list.
  if (offsets.size() - num_offset == 1)
  {
    offsets[num_offset++] = num_adj;
  }
  assert( (unsigned)num_offset == offsets.size() );
  
  if (sizeof(size_t) != sizeof(int))  // do copy if could not be avoided
    copy_from_sidl( sidl_offsets, &offsets[0] );
}


void MeshTSTTImpl::vertices_get_attached_elements( 
                                     const VertexHandle* vertices,
                                     size_t num_vertex,
                                     msq_std::vector<ElementHandle>& elements,
                                     msq_std::vector<size_t>& offsets,
                                     MsqError& err )
{
  get_adjacent_entities( vertices, num_vertex, inputSetType, elements, offsets, err ); 
  MSQ_ERRRTN(err);
  
    // We need to use inputElements.  Fill it if it hasn't been filled yet.
  if (inputElements.empty())
  {
    populate_input_elements(err);
    MSQ_ERRRTN(err);
  }
  
    // Remove all elements not in inputElements
  msq_std::vector<size_t>::iterator offset_iter = offsets.begin();
  size_t read_idx, write_idx;
  for (read_idx = write_idx = 0; read_idx < elements.size(); ++read_idx)
  {
    if (*offset_iter == read_idx)
    {
      *offset_iter = write_idx;
      ++offset_iter;
    }

    if (msq_std::binary_search( inputElements.begin(), inputElements.end(), elements[read_idx] ))
      elements[write_idx++] = elements[read_idx];
  }
  assert( offset_iter + 1 == offsets.end() && *offset_iter == read_idx );
  *offset_iter = write_idx;
}


//**************** Element Topology *****************


/** Get connectivity
 *\param elements - Array of length num_elems containing elements
 *                  handles of elements for which connectivity is to
 *                  be queried.
 *\param vertices - Array of vertex handles in connectivity list.
 *\param offsets  - Indices into \ref vertex_handles, one per element
 */
void MeshTSTTImpl::elements_get_attached_vertices(
  const ElementHandle *elements,
  size_t num_elems,
  msq_std::vector<VertexHandle>& vertices,
  msq_std::vector<size_t>& offsets,
  Mesquite::MsqError &err)
{
  get_adjacent_entities( elements, num_elems, TSTTM::EntityType_VERTEX, vertices, offsets, err );
  MSQ_CHKERR(err);
}


void MeshTSTTImpl::get_all_elements( msq_std::vector<ElementHandle>& elements,
                                     MsqError& err )
{
  try {
    int count = meshIFace.getNumOfType( elementSet, TSTTM::EntityType_ALL_TYPES );
    elements.resize( count );
    sidl::array<void*> handles( convert_to_sidl_vector( &elements[0], count ) );
    meshIFace.getEntities( elementSet, 
                           TSTTM::EntityType_ALL_TYPES, 
                           TSTTM::EntityTopology_ALL_TOPOLOGIES,
                           handles, count );
    assert( (unsigned)count == elements.size() );
  } 
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
}

void MeshTSTTImpl::get_all_vertices( msq_std::vector<VertexHandle>& vertices,
                                     MsqError& err )
{
  try {
    int count = meshIFace.getNumOfType( nodeSet, TSTTM::EntityType_ALL_TYPES );
    vertices.resize( count );
    sidl::array<void*> handles( convert_to_sidl_vector( &vertices[0], count ) );
    meshIFace.getEntities( nodeSet, 
                           TSTTM::EntityType_ALL_TYPES, 
                           TSTTM::EntityTopology_ALL_TOPOLOGIES,
                           handles, count );
    assert( (unsigned)count == vertices.size() );
  } 
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
}
      

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
void MeshTSTTImpl::elements_get_topologies(
  const ElementHandle *element_handle_array,
  EntityTopology *element_topologies,
  size_t num_elements, MsqError &err)
{
  try {
    
    sidl::array<TSTTM::EntityTopology> topologies( alloc_sidl_vector<TSTTM::EntityTopology>( num_elements ) );
    sidl::array<void*> handles( convert_to_sidl_vector( const_cast<void**>(element_handle_array), num_elements ) );
    int topo_len_out;
    arrIFace.getEntArrTopo( handles, num_elements, topologies, topo_len_out );
    if ((size_t)topo_len_out != num_elements) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return;
    }
    
    for (unsigned i = 0; i < num_elements; ++i)
      element_topologies[i] = topologyMap[topologies.get(i)];
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void MeshTSTTImpl::release_entity_handles(
  const Mesquite::Mesh::EntityHandle */*handle_array*/,
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
  TSTTB::TagValueType tstt_type;
  switch (type) {
    case Mesquite::Mesh::BYTE:   tstt_type = TSTTB::TagValueType_BYTES;         break;
    case Mesquite::Mesh::INT:    tstt_type = TSTTB::TagValueType_INTEGER;       break;
    case Mesquite::Mesh::DOUBLE: tstt_type = TSTTB::TagValueType_DOUBLE;        break;
    case Mesquite::Mesh::HANDLE: tstt_type = TSTTB::TagValueType_ENTITY_HANDLE; break;
    default:
      MSQ_SETERR(err)("Invalid tag type", MsqError::INVALID_ARG );
      return 0;
  }
  
  try {
    void* handle = 0;
    tagIFace.createTag( name, length, tstt_type, handle );
    return handle;
  } 
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

void MeshTSTTImpl::tag_delete( TagHandle handle, MsqError& err )
{
  try {
    tagIFace.destroyTag( handle, true );
  } 
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

TagHandle MeshTSTTImpl::tag_get( const msq_std::string& name, MsqError& err )
{
  try {
    return tagIFace.getTagHandle( name );
  } 
  catch(::TSTTB::Error &tstt_err) {
    if (tstt_err.getErrorType() != TSTTB::ErrorType_TAG_NOT_FOUND) {
      MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    } 
    else {
      MSQ_SETERR(err)( name, MsqError::TAG_NOT_FOUND );
    }
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
  TSTTB::TagValueType type;
  try {
    name = tagIFace.getTagName( handle );
    size = tagIFace.getTagSizeBytes( handle );
    type = tagIFace.getTagType( handle );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
  
  size_t tsize;
  switch (type) {
    case TSTTB::TagValueType_BYTES        : tsize = sizeof(char  ); type_out = Mesquite::Mesh::BYTE  ; break;
    case TSTTB::TagValueType_INTEGER      : tsize = sizeof(int   ); type_out = Mesquite::Mesh::INT   ; break;
    case TSTTB::TagValueType_DOUBLE       : tsize = sizeof(double); type_out = Mesquite::Mesh::DOUBLE; break;
    case TSTTB::TagValueType_ENTITY_HANDLE: tsize = sizeof(void* ); type_out = Mesquite::Mesh::HANDLE; break;
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
    size_t len, size = tagIFace.getTagSizeBytes( tag );
    int count = 0;
    sidl::array<void*> handles( convert_to_sidl_vector( const_cast<void**>(array), num_elems ) );
    switch (tagIFace.getTagType( tag ))
    {
      case TSTTB::TagValueType_ENTITY_HANDLE:
      {
        len = size / sizeof(void*);
        sidl::array<void*> sdata( convert_to_sidl_vector( (void**)data, len*num_elems ));
        arrTagIFace.setEHArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_DOUBLE:
      {
        len = size / sizeof(double);
        sidl::array<double> sdata( convert_to_sidl_vector( (double*)data, len*num_elems ));
        arrTagIFace.setDblArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_INTEGER:
      {
        len = size / sizeof(int);
        sidl::array<int> sdata( convert_to_sidl_vector( (int*)data, len*num_elems ));
        arrTagIFace.setIntArrData( handles, num_elems, tag, sdata, count );
      }
      break;

      default:
      {
        len = size * num_elems;
        sidl::array<char> sdata( convert_to_sidl_vector( (char*)data, len ) );
        arrTagIFace.setArrData( handles, num_elems, tag, sdata, count );
      }
    }
  }
  catch(::TSTTB::Error &tstt_err) {
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
    int count;
    size_t len, size = tagIFace.getTagSizeBytes( tag );
    sidl::array<void*> handles = convert_to_sidl_vector( (void**)array, num_elems );
    switch (tagIFace.getTagType( tag ))
    {
      case TSTTB::TagValueType_ENTITY_HANDLE:
      {
        len = size / sizeof(void*);
        sidl::array<void*> sdata( convert_to_sidl_vector( (void**)data, len*num_elems ));
        arrTagIFace.getEHArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_DOUBLE:
      {
        len = size / sizeof(double);
        sidl::array<double> sdata( convert_to_sidl_vector( (double*)data, len*num_elems ));
        arrTagIFace.getDblArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_INTEGER:
      {
        len = size / sizeof(int);
        sidl::array<int> sdata( convert_to_sidl_vector( (int*)data, len*num_elems ));
        arrTagIFace.getIntArrData( handles, num_elems, tag, sdata, count );
      }
      break;

      default:
      {
        len = size * num_elems;
        sidl::array<char> sdata = convert_to_sidl_vector( (char*)data, len );
        arrTagIFace.getArrData( handles, num_elems, tag, sdata, count );
      }
    }
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

} // namespace Mesquite
