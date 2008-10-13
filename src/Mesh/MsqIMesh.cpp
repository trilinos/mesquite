/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MsqIMesh.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "MsqIMesh.hpp"
#include "MsqError.hpp"
#include "MeshInterface.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include <assert.h>
#include "MsqIBase.hpp"
#include <algorithm>

namespace Mesquite {

class iMeshArrIter : public EntityIterator
{
  private:
  
    std::vector<iBase_EntityHandle> handleArray;
    iMesh_EntityArrIterator imeshIter;
    int index, count, notAtEnd;
    iMesh_Instance theMesh;
    
    inline void get_next_array( int& err )
    {
      int alloc = handleArray.size();
      index = count = 0;
      iBase_EntityHandle* ptr = &handleArray[0];
      iMesh_getNextEntArrIter( theMesh, imeshIter, &ptr, &alloc, &count, &notAtEnd, &err );
    }
    
  public:
  
    iMeshArrIter( iMesh_Instance mesh,
                  iBase_EntitySetHandle meshset,
                  int type, int topo, int& err,
                  unsigned buffer_count = 1024 );
    
    virtual void restart();
    virtual iBase_EntityHandle operator*() const;
    virtual bool is_at_end() const;
    virtual void operator++();
    
    virtual ~iMeshArrIter();
};

iMeshArrIter::iMeshArrIter( iMesh_Instance mesh,
                            iBase_EntitySetHandle meshset,
                            int type, int topo, int& err,
                            unsigned buffer_count )
  : handleArray(buffer_count),
    imeshIter(0),
    index(0), count(0),
    notAtEnd(0),
    theMesh(mesh)
{
  iMesh_initEntArrIter( theMesh, meshset, type, topo, buffer_count, &imeshIter, &err );
  if (err == iBase_SUCCESS)
    get_next_array( err );
  else
    imeshIter = 0;
}

iMeshArrIter::~iMeshArrIter()
{
  int err;
  if (imeshIter)
    iMesh_endEntArrIter( theMesh, imeshIter, &err );
}

void iMeshArrIter::restart()
{
  int err;
  iMesh_resetEntArrIter( theMesh, imeshIter, &err );
}

iBase_EntityHandle iMeshArrIter::operator*() const
{
  return handleArray[index];
}

bool iMeshArrIter::is_at_end() const
{
  return index == count && !notAtEnd;
}

void iMeshArrIter::operator++() 
{
  ++index;
  if (index == count && notAtEnd) {
    int err;
    get_next_array( err );
  }
}


class MsqIMeshImpl : public MsqIMesh
{
  public:

    MsqIMeshImpl(iMesh_Instance mesh, 
                 const char* fixed_tag_name,
                 MsqError &err);
    virtual ~MsqIMeshImpl();
    
      /** \brief set mesh to be smoothed.
       *
       * Set the mesh which Mesquite is to smooth.  Optionally
       * specify fixed vertices.
       * NOTE: If an active set is not specified, the default
       *       is to use the global set (the ENTIRE mesh.)
       *
       *\param element_set ITAPS entity set handle for set containing
       *                  mesh elements and vertices for which quality 
       *                  is to be improved.
       */
    virtual void set_active_set( void* element_set, 
                                 iBase_EntityType type,
                                 MsqError& );
  
    virtual iMesh_Instance get_imesh_instance() const;
    virtual iBase_EntitySetHandle get_entity_set() const;
    

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
                                iBase_EntityType target_type,
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
    
    /** The IMesh instance */
    iMesh_Instance meshInstance;
    
    /** Have mesh */
    bool haveMesh;
    /** ITAPS entity set handle for elements to improve */
    iBase_EntitySetHandle elementSet;
    /** ITAPS entity set handle for nodes to move */
    iBase_EntitySetHandle nodeSet;
    /** std::set containing elements in elementSet, used
     *  to constrain vertex->element adjaceny queries to
     *  only those elements that are in the input element set.
     */
    msq_std::vector<iBase_EntityHandle> inputElements;
    
    /** The type of elements contained in the input element set.
     * Should be one of:
     * - iBase_REGION    - volume elements
     * - iBase_FACE      - face/2d elements
     * - iBase_ALL_TYPES - mixed volume and face elements
     */
    iBase_EntityType inputSetType;
    
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
    /** Dimension is queried once during create and cached */
    int geometricDimension;
    /** Map iMesh_EntityTopology to Mesquite::EntityTopology */
    EntityTopology topologyMap[iMesh_ALL_TOPOLOGIES+1];
};

/*************************************************************************
 *                          Mesh Definition
 ************************************************************************/

MsqIMesh* MsqIMesh::create( iMesh_Instance mesh, 
                            iBase_EntitySetHandle meshset, 
                            iBase_EntityType type,
                            MsqError& err,
                            const char* fixed_tag_name )
{
  MsqIMesh* result = new MsqIMeshImpl( mesh, fixed_tag_name, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  result->set_active_set( meshset, type, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  return result;
}

MsqIMesh* MsqIMesh::create( iMesh_Instance mesh, 
                            MsqError& err,
                            const char* fixed_tag_name )
{
  MsqIMesh* result = new MsqIMeshImpl( mesh, fixed_tag_name, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  return result;
}

MsqIMesh::~MsqIMesh() {}


MsqIMeshImpl::MsqIMeshImpl( iMesh_Instance itaps_mesh, 
                            const char* fixed_tag_name,
                            Mesquite::MsqError& err ) 
  : meshInstance(itaps_mesh), 
    elementSet(0), nodeSet(0), 
    inputSetType( iBase_ALL_TYPES ),
    byteTag(0), createdByteTag(false),
    fixedTag(0), createdFixedTag(false),
    geometricDimension(0)
{
    // Initialize topology map 
  
  const size_t mapsize = sizeof(topologyMap) / sizeof(Mesquite::EntityTopology);
  if (mapsize < iMesh_ALL_TOPOLOGIES)
  {
    MSQ_SETERR(err)("MsqIMesh needs to be updated for new iMesh element topologies.",
                    MsqError::INTERNAL_ERROR);
    return;
  }
  
  for (size_t i = 0; i <= iMesh_ALL_TOPOLOGIES; ++i)
    topologyMap[i] = Mesquite::MIXED;
  
  topologyMap[iMesh_TRIANGLE     ] = Mesquite::TRIANGLE;
  topologyMap[iMesh_QUADRILATERAL] = Mesquite::QUADRILATERAL;
  topologyMap[iMesh_TETRAHEDRON  ] = Mesquite::TETRAHEDRON;
  topologyMap[iMesh_HEXAHEDRON   ] = Mesquite::HEXAHEDRON;
  topologyMap[iMesh_PRISM        ] = Mesquite::PRISM;
  topologyMap[iMesh_PYRAMID      ] = Mesquite::PYRAMID;
  
      // Get tag for fixed flag
  if (fixed_tag_name == 0)
    fixed_tag_name = VERTEX_FIXED_TAG_NAME;
  int ierr;
  iMesh_getTagHandle( meshInstance, fixed_tag_name, 
                      &fixedTag, &ierr,
                      strlen(fixed_tag_name) );
  if (iBase_SUCCESS == ierr) {
    int size, type;
    iMesh_getTagSizeBytes( meshInstance, fixedTag, &size, &ierr );
    if (iBase_SUCCESS != ierr || size != sizeof(int)) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE,
                       "Tag \"%s\" exists with invalid size", 
                       fixed_tag_name );
      return;
    }
    iMesh_getTagType( meshInstance, fixedTag, &type, &ierr );
    if (iBase_SUCCESS != ierr || type != iBase_INTEGER) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE,
                       "Tag \"%s\" exists with invalid type", 
                       fixed_tag_name );
      return;
    }
  }
  else {
    fixedTag = 0;
  }
  
    // Get/create tag for vertex byte
  iMesh_getTagHandle( meshInstance, 
                      VERTEX_BYTE_TAG_NAME,
                      &byteTag, &ierr,
                      strlen(VERTEX_BYTE_TAG_NAME) );
  if (iBase_SUCCESS != ierr) {
    iMesh_createTag( meshInstance, 
                     VERTEX_BYTE_TAG_NAME,
                     1, iBase_INTEGER,
                     &byteTag, &ierr,
                     strlen(VERTEX_BYTE_TAG_NAME) );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                       "Tag \"%s\" could not be created", 
                       VERTEX_BYTE_TAG_NAME );
      return;
    }
  }
  else {
    int size, type;
    iMesh_getTagSizeBytes( meshInstance, byteTag, &size, &ierr );
    if (iBase_SUCCESS != ierr || size != sizeof(int)) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE,
                       "Tag \"%s\" exists with invalid size", 
                       VERTEX_BYTE_TAG_NAME );
      return;
    }
    iMesh_getTagType( meshInstance, byteTag, &type, &ierr );
    if (iBase_SUCCESS != ierr || type != iBase_INTEGER) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE,
                       "Tag \"%s\" exists with invalid type", 
                       VERTEX_BYTE_TAG_NAME );
      return;
    }
  }
  iMesh_getGeometricDimension( meshInstance, &geometricDimension, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}

MsqIMeshImpl::~MsqIMeshImpl() 
{
  int err;
  if (elementSet) {
    iMesh_destroyEntSet( meshInstance, elementSet, &err );
    if (iBase_SUCCESS != err)
      process_itaps_error(err);
  }
  if (nodeSet) {
    iMesh_destroyEntSet( meshInstance, nodeSet, &err );
    if (iBase_SUCCESS != err)
      process_itaps_error(err);
  }
}


void MsqIMeshImpl::set_int_tag( iBase_TagHandle tag,
                                iBase_EntitySetHandle elem_set, 
                                int value, 
                                MsqError& err )
{
  const unsigned BUFFER_COUNT = 1024;
  iBase_EntityHandle handle_array[BUFFER_COUNT];
  int value_array[BUFFER_COUNT];
  if (!value)
    memset( value_array, 0, sizeof(value_array ));
  else for (unsigned i = 0; i < BUFFER_COUNT; ++i)
    value_array[i] = value;
  
  iMesh_EntityArrIterator iter;
  int more = 1, ierr, alloc = BUFFER_COUNT;
  
  assert(elem_set);
  iMesh_initEntArrIter( meshInstance, elem_set, iBase_VERTEX, iMesh_POINT, BUFFER_COUNT, &iter, &ierr );
  if (ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  for (;;) {
    int count = 0;
    iBase_EntityHandle* ptr = handle_array;
    iMesh_getNextEntArrIter( meshInstance, iter, &ptr, &alloc, &count, &more, &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      iMesh_endEntArrIter( meshInstance, iter, &ierr );
      return;
    }
    if (!count)
      break;
    iMesh_setIntArrData( meshInstance, handle_array, count, tag, value_array, count, &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      iMesh_endEntArrIter( meshInstance, iter, &ierr );
      return;
    }
  }
 
  iMesh_endEntArrIter( meshInstance, iter, &ierr );
  if (ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}

iMesh_Instance MsqIMeshImpl::get_imesh_instance() const
{
  return meshInstance;
}

iBase_EntitySetHandle MsqIMeshImpl::get_entity_set() const
{
  return elementSet;
}

void MsqIMeshImpl::set_active_set( iBase_EntitySetHandle elem_set, 
                                   iBase_EntityType type_in,
                                   MsqError& err )
{
  const int ELEM_BUFFER_SIZE = 1024;
  const int NODE_BUFFER_SIZE = 27 * ELEM_BUFFER_SIZE; 
  iBase_EntityHandle elements[ELEM_BUFFER_SIZE], nodes[NODE_BUFFER_SIZE];
  int offsets[ELEM_BUFFER_SIZE+1], ierr;
  iMesh_EntityArrIterator iter = 0;
 
  if (elementSet)
    iMesh_destroyEntSet( meshInstance, elementSet, &ierr );
    
  if (nodeSet)
    iMesh_destroyEntSet( meshInstance, nodeSet, &ierr );
  
  iMesh_createEntSet( meshInstance, 0, &elementSet, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  iMesh_createEntSet( meshInstance, 0, &nodeSet, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
    
    // Iterate over set twice, once for FACEs and once for REGIONs
  bool have_faces = false, have_regions = false;
  int start = (type_in == iBase_ALL_TYPES) ? 0 : 1; // don't loop twice if user specified type
  for (int i = start; i < 2; ++i)
  {
    iBase_EntityType type = type_in;
    if (type == iBase_ALL_TYPES)
      type = i ? iBase_REGION : iBase_FACE;
    bool& have_some = (type == iBase_REGION) ? have_regions : have_faces;
    
    iMesh_initEntArrIter( meshInstance, elem_set, type, iMesh_ALL_TOPOLOGIES, ELEM_BUFFER_SIZE, &iter, &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      return;
    }
    

    for (;;) {
      int count = 0;
      int junk = ELEM_BUFFER_SIZE, junk2;
      iBase_EntityHandle* ptr3 = elements;
      iMesh_getNextEntArrIter( meshInstance, iter, &ptr3, &junk, &count, &junk2, &ierr );
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
        return;
      }

      if (!count)
        break;
      have_some = true;
      
      iMesh_addEntArrToSet( meshInstance, elements, count, &elementSet, &ierr );
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
        return;
      }
     
      int num_nodes, junk3;
      junk = NODE_BUFFER_SIZE;
      junk2 = ELEM_BUFFER_SIZE + 1;
      int* ptr = offsets;
      iBase_EntityHandle* ptr2 = nodes;
      iMesh_getEntArrAdj( meshInstance,
                          elements, count,
                          iBase_VERTEX,
                          &ptr2, &junk, &num_nodes,
                          &ptr, &junk2, &junk3, 
                          &ierr );
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
        return;
      }

      iMesh_addEntArrToSet( meshInstance, nodes, num_nodes, &nodeSet, &ierr );
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
        return;
      }
    }
    
    iMesh_endEntArrIter( meshInstance, iter, &ierr );
    iter = 0;
  } // for (type)

  if (!have_faces)
    inputSetType = iBase_REGION;
  else if (!have_regions)
    inputSetType = iBase_FACE;
  else
    inputSetType = iBase_ALL_TYPES;
  
    // clear cached data
  inputElements.clear();

    // Clear vertex byte tag
  set_int_tag( byteTag, nodeSet, 0, err );

  MSQ_CHKERR(err);
}

void MsqIMeshImpl::populate_input_elements( MsqError& err ) 
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
int MsqIMeshImpl::get_geometric_dimension(Mesquite::MsqError &err)
{
  return geometricDimension;
}
    
    
// Returns a pointer to an iterator that iterates over the
// set of all vertices in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If vertices are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
VertexIterator* MsqIMeshImpl::vertex_iterator(MsqError& err)
{
  int ierr;
  VertexIterator* iter = new iMeshArrIter( meshInstance, 
                                           nodeSet, 
                                           iBase_ALL_TYPES,
                                           iMesh_ALL_TOPOLOGIES,
                                           ierr );
  if (iBase_SUCCESS != ierr) {
    delete iter;
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return 0;
  }
  
  return iter;
}
    
// Returns a pointer to an iterator that iterates over the
// set of all top-level elements in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If elements are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
ElementIterator* MsqIMeshImpl::element_iterator(MsqError &err)
{
  int ierr;
  VertexIterator* iter = new iMeshArrIter( meshInstance, 
                                           elementSet, 
                                           iBase_ALL_TYPES,
                                           iMesh_ALL_TOPOLOGIES,
                                           ierr );
  if (iBase_SUCCESS != ierr) {
    delete iter;
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return 0;
  }
  
  return iter;
}

//************ Vertex Properties ********************
// Returns true or false, indicating whether the vertex
// is allowed to be repositioned.  True indicates that the vertex
// is fixed and cannot be moved.  Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
void MsqIMeshImpl::vertices_get_fixed_flag(
  const VertexHandle vert_array[], 
  bool bool_array[],
  size_t num_vtx, MsqError &err)
{
    // If mesh does not contain a fixed tag, assume no vertices are fixed
  if (!fixedTag) {
    memset( bool_array, 0, num_vtx * sizeof(bool) );
    return;
  }
  
  std::vector<int> values(num_vtx);
  int ierr, junk = num_vtx, junk2 = num_vtx;
  int* ptr = &values[0];
  iMesh_getIntArrData( meshInstance, vert_array, num_vtx, fixedTag, &ptr, &junk, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  for (size_t i = 0; i < num_vtx; ++i)
    bool_array[i] = values[i];
}

// Get vertex coordinates 
void MsqIMeshImpl::vertices_get_coordinates(
  const Mesquite::Mesh::VertexHandle vert_array[],
  MsqVertex* coordinates, 
  size_t num_vtx, 
  MsqError &err)
{
  int order = iBase_UNDETERMINED;
  std::vector<double> dbl_store( 3*num_vtx );
  double* dbl_array = &dbl_store[0];
  
  int ierr, junk = 3*num_vtx, junk2;
  iMesh_getVtxArrCoords( meshInstance, vert_array, num_vtx, &order, &dbl_array, &junk, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  if (geometricDimension == 2)
  {
    if (order == iBase_INTERLEAVED)
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
  else 
  {
    if (order == iBase_INTERLEAVED)
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
}

void MsqIMeshImpl::vertex_set_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  const Vector3D &coords, MsqError &err)
{
  int ierr;
  iMesh_setVtxCoords( meshInstance, vertex, coords[0], coords[1], coords[2], &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
}

// Each vertex has a byte-sized flag that can be used to store
// flags.  This byte's value is neither set nor used by the mesh
// implementation.  It is intended to be used by Mesquite algorithms.
// Until a vertex's byte has been explicitly set, its value is 0.
void MsqIMeshImpl::vertex_set_byte (
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char byte, MsqError &err)
{
  int ierr, value = byte;
  iMesh_setIntData( meshInstance, vertex, byteTag, value, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
}

void MsqIMeshImpl::vertices_set_byte (
  const VertexHandle *vert_array,
  const unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
  std::vector<int> data(array_size);
  std::copy( byte_array, byte_array + array_size, data.begin() );
  int ierr;
  iMesh_setIntArrData( meshInstance, const_cast<VertexHandle*>(vert_array), array_size, byteTag, &data[0], array_size, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
}

// Retrieve the byte value for the specified vertex or vertices.
// The byte value is 0 if it has not yet been set via one of the
// *_set_byte() functions.
void MsqIMeshImpl::vertex_get_byte(
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char *byte, MsqError &err)
{
  int ierr, value;
  iMesh_getIntData( meshInstance, vertex, byteTag, &value, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
  *byte = value;
}

void MsqIMeshImpl::vertices_get_byte(
  const VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
  std::vector<int> data(array_size);
  int ierr;
  int* ptr = &data[0];
  int junk1 = data.size(), junk2;
  iMesh_getIntArrData( meshInstance, const_cast<VertexHandle*>(vert_array), array_size, byteTag, &ptr, &junk1, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
  std::copy( data.begin(), data.end(), byte_array );
}


//**************** Topology *****************

void MsqIMeshImpl::get_adjacent_entities( const iBase_EntityHandle* source,
                                          size_t num_source,
                                          iBase_EntityType target_type,
                                          msq_std::vector<iBase_EntityHandle>& target,
                                          msq_std::vector<size_t>& offsets,
                                          MsqError& err )
{
  if (num_source == 0)
    return;
  
  int ierr, num_adj = 0, num_offset;
  
  assert( sizeof(size_t) >= sizeof(int) );
  offsets.resize( num_source + 1 );
  int* ptr2 = (int*)&offsets[0];
  bool expand = false;
  if (sizeof(size_t) > sizeof(int))
    expand = true;
  
  bool have_adj = false;
    // If passed vector has allocated storage, try to use existing space
  if (target.capacity() >= num_source)
  {
    target.resize( target.capacity() );
    int junk1 = target.capacity(), junk3 = offsets.size();
    iBase_EntityHandle* ptr = &target[0];
    iMesh_getEntArrAdj( meshInstance, source, num_source,
                        target_type, 
                        &ptr, &junk1, &num_adj, 
                        &ptr2, &junk3, &num_offset, 
                        &ierr );
    if (iBase_SUCCESS == ierr) {
      have_adj = true;
      target.resize( num_adj );
    }
  }
  
    // If implementation passed back a size, try that
  if (!have_adj && num_adj && (unsigned)num_adj > target.capacity())
  {
    target.resize( num_adj );
    int junk1 = target.capacity(), junk3 = offsets.size();
    iBase_EntityHandle* ptr = &target[0];
    iMesh_getEntArrAdj( meshInstance, source, num_source,
                        target_type, 
                        &ptr, &junk1, &num_adj, 
                        &ptr2, &junk3, &num_offset, 
                        &ierr );
    if (iBase_SUCCESS == ierr)
      have_adj = true;
  }

    // Try with empty sidl array, and copy into elements vector
  if (!have_adj)
  {
    iBase_EntityHandle* mArray = 0;
    int junk1 = 0, junk3 = offsets.size();
    iMesh_getEntArrAdj( meshInstance, source, num_source,
                        target_type, 
                        &mArray, &junk1, &num_adj, 
                        &ptr2, &junk3, &num_offset, 
                        &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      return;
    }
    
    target.resize( num_adj );
    std::copy( mArray, mArray + num_adj, target.begin() );
    free( mArray );
  }
  
  if (expand) {
    for (size_t i = num_offset; i > 0; --i)
      offsets[i-1] = ptr2[i-1];
  }
  
  // iMesh implementations seem to be inconsistent with regard to 
  // placing the last value on this list.
  if (offsets.size() - num_offset == 1)
    offsets[num_offset++] = num_adj;
  assert( (unsigned)num_offset == offsets.size() );
}


void MsqIMeshImpl::vertices_get_attached_elements( 
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
void MsqIMeshImpl::elements_get_attached_vertices(
  const ElementHandle *elements,
  size_t num_elems,
  msq_std::vector<VertexHandle>& vertices,
  msq_std::vector<size_t>& offsets,
  Mesquite::MsqError &err)
{
  get_adjacent_entities( elements, num_elems, iBase_VERTEX, vertices, offsets, err );
  MSQ_CHKERR(err);
}


void MsqIMeshImpl::get_all_elements( msq_std::vector<ElementHandle>& elements,
                                     MsqError& err )
{
  int count, ierr;
  iMesh_getNumOfType( meshInstance, elementSet, iBase_ALL_TYPES, &count, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  elements.resize(count);
  int junk = count;
  iBase_EntityHandle* ptr = &elements[0];
  iMesh_getEntities( meshInstance,
                     elementSet,
                     iBase_ALL_TYPES,
                     iMesh_ALL_TOPOLOGIES,
                     &ptr, &junk, &count,
                     &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  assert( (unsigned)count == elements.size() );
}

void MsqIMeshImpl::get_all_vertices( msq_std::vector<VertexHandle>& vertices,
                                     MsqError& err )
{
  int count, ierr;
  iMesh_getNumOfType( meshInstance, nodeSet, iBase_ALL_TYPES, &count, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  vertices.resize(count);
  int junk = count;
  iBase_EntityHandle* ptr = &vertices[0];
  iMesh_getEntities( meshInstance,
                     nodeSet,
                     iBase_ALL_TYPES,
                     iMesh_ALL_TOPOLOGIES,
                     &ptr, &junk, &count,
                     &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  assert( (unsigned)count == vertices.size() );
}
      

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
void MsqIMeshImpl::elements_get_topologies(
  const ElementHandle *element_handle_array,
  EntityTopology *element_topologies,
  size_t num_elements, MsqError &err)
{
    // don't copy unless we have to
  std::vector<int> topo_store;
  int* topo_array;
  if (sizeof(EntityTopology) == sizeof(int))
    topo_array = (int*)element_topologies;
  else {
    topo_store.resize(num_elements);
    topo_array = &topo_store[0];
  }
  
  int ierr, junk1 = num_elements, junk2;
  iMesh_getEntArrTopo( meshInstance, element_handle_array, num_elements, &topo_array, &junk1, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  for (size_t i = 0; i < num_elements; ++i)
    element_topologies[i] = topologyMap[topo_array[i]];
}

//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void MsqIMeshImpl::release_entity_handles(
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
void MsqIMeshImpl::release()
{
}

//**************** Tags ****************
TagHandle MsqIMeshImpl::tag_create( const msq_std::string& name, 
                                    TagType type, unsigned length,
                                    const void* ,
                                    MsqError& err )
{
  int itaps_type;
  switch (type) {
    case Mesquite::Mesh::BYTE:   itaps_type = iBase_BYTES;         break;
    case Mesquite::Mesh::INT:    itaps_type = iBase_INTEGER;       break;
    case Mesquite::Mesh::DOUBLE: itaps_type = iBase_DOUBLE;        break;
    case Mesquite::Mesh::HANDLE: itaps_type = iBase_ENTITY_HANDLE; break;
    default:
      MSQ_SETERR(err)("Invalid tag type", MsqError::INVALID_ARG );
      return 0;
  }
  
  int ierr;
  TagHandle result;
  iMesh_createTag( meshInstance, 
                   name.c_str(), 
                   length, itaps_type,
                   &result, &ierr,
                   name.size() );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return 0;
  }
  
  return result;
}

void MsqIMeshImpl::tag_delete( TagHandle handle, MsqError& err )
{
  int ierr;
  iMesh_destroyTag( meshInstance, handle, true, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}

TagHandle MsqIMeshImpl::tag_get( const msq_std::string& name, MsqError& err )
{
  TagHandle handle = 0;
  int ierr;
  iMesh_getTagHandle( meshInstance, name.c_str(), &handle, &ierr, name.length() );
  if (iBase_TAG_NOT_FOUND == ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::TAG_NOT_FOUND );
    return 0;
  }
  else if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return 0;
  }
  return handle;
}


void MsqIMeshImpl::tag_properties( TagHandle handle,
                                   msq_std::string& name_out,
                                   TagType& type_out,
                                   unsigned& length_out,
                                   MsqError& err )
{
  char buffer[256];
  int ierr1, ierr2, ierr3, itype;
  
  iMesh_getTagName( meshInstance, handle, buffer, &ierr1, sizeof(buffer) );
  iMesh_getTagSizeValues( meshInstance, handle, (int*)&length_out, &ierr2 );
  iMesh_getTagType( meshInstance, handle, &itype, &ierr3 );
  
  int ierr = iBase_SUCCESS;
  if (ierr1 != iBase_SUCCESS)
    ierr = ierr1;
  else if (ierr2 != iBase_SUCCESS)
    ierr = ierr2;
  else if (ierr3 != iBase_SUCCESS)
    ierr = ierr3;
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  buffer[255] = '\0';
  name_out = buffer;
  switch (itype) {
    case iBase_BYTES        : type_out = Mesquite::Mesh::BYTE  ; break;
    case iBase_INTEGER      : type_out = Mesquite::Mesh::INT   ; break;
    case iBase_DOUBLE       : type_out = Mesquite::Mesh::DOUBLE; break;
    case iBase_ENTITY_HANDLE: type_out = Mesquite::Mesh::HANDLE; break;
    default:
      MSQ_SETERR(err)("Unsupported iMesh tag type", MsqError::NOT_IMPLEMENTED );
      return;
  }
}

void MsqIMeshImpl::tag_set_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         const void* data,
                                         MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}

void MsqIMeshImpl::tag_set_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        const void* data,
                                        MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}
    
void MsqIMeshImpl::tag_set_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 const void* data,
                                 MsqError& err )
{
  int ierr, size;
  iMesh_getTagSizeBytes( meshInstance, tag, &size, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  iMesh_setArrData( meshInstance, const_cast<iBase_EntityHandle*>(array), num_elems, tag, static_cast<const char*>(data), size*num_elems, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}  


void MsqIMeshImpl::tag_get_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         void* data,
                                         MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}

void MsqIMeshImpl::tag_get_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        void* data,
                                        MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}
    
void MsqIMeshImpl::tag_get_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 void* data,
                                 MsqError& err )
{
  int ierr, size;
  iMesh_getTagSizeBytes( meshInstance, tag, &size, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  char* ptr = static_cast<char*>(data);
  int junk1 = size*num_elems, junk2;
  iMesh_getArrData( meshInstance, array, num_elems, tag, &ptr, &junk1, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}


} // namespace Mesquite
