/*!
  \file   MeshTSTT.cpp
  \brief  


  \author Thomas Leurent
  \date   2003-06-12
*/

#include <set>
//#include <pair>

#include "MeshTSTT.hpp"
#include "TSTT.hh"

using std::cout;
using std::cerr;
using std::endl;

namespace
{

#define ENTIRE_MESH 0  // passing the null pointer as the first argument of
                   // TSTT::entitySet functions queries the data for the whole mesh
  
  class MeshTSTT_EntityIterator : public Mesquite::EntityIterator
  {
  public:
    MeshTSTT_EntityIterator(TSTT::LocalTSTTMesh& tstt_mesh,
                            const TSTT::EntityType& entity_type)
      : iteratorEntityType(entity_type)
    {
      iteratorMesh = tstt_mesh;
      iteratorMesh.entitysetInitializeWorksetIterator(
                                         ENTIRE_MESH,
                                         iteratorEntityType,
                                         ::TSTT::ALL_TOPOLOGIES,
                                         1, // requested workset size
                                         tsttWorksetIterator);
      entityHandle = ::SIDL::array<void*>::create1d(1); // array with one entry.
      bool more = iteratorMesh.entitysetGetNextWorkset(tsttWorksetIterator, entityHandle);
      if (!more)
        entityHandle.set(0, NULL); 
    }
    
      // Moves the iterator back to the first
      // entity in the list.
    virtual void restart()
      {
        iteratorMesh.entitysetDestroyWorksetIterator(
                                 tsttWorksetIterator);
        iteratorMesh.entitysetInitializeWorksetIterator(
                                         ENTIRE_MESH,
                                         iteratorEntityType,
                                         ::TSTT::ALL_TOPOLOGIES,
                                         1, // requested workset size
                                         tsttWorksetIterator);
      }
    
      // *iterator.  Return the handle currently
      // being pointed at by the iterator.
    virtual Mesquite::Mesh::EntityHandle operator*() const
      {
        return entityHandle.get(0); // NULL if past the end. 
      }
    
      // ++iterator
    virtual void operator++()
      {
        bool more = iteratorMesh.entitysetGetNextWorkset(tsttWorksetIterator, entityHandle);
        if (!more)
          entityHandle.set(0, NULL); 
      }
    
      // Returns false until the iterator has
      // been advanced PAST the last entity.
      // Once is_at_end() returns true, *iterator
      // returns NULL.
    virtual bool is_at_end() const
      {
        return ( entityHandle.get(0)!=NULL ? false : true);
      }
    
  private:
    void* tsttWorksetIterator;
    ::TSTT::LocalTSTTMesh iteratorMesh;
    ::TSTT::EntityType iteratorEntityType;
    ::SIDL::array<void*> entityHandle; // Array has one entry only. That entry is
                                     // set to NULL if the end of the mesh is reached. 
  };

  
#undef __FUNC__
#define __FUNC__ "::mesquite_equivalent_topology"
  inline Mesquite::EntityTopology mesquite_equivalent_topology(
                                     const TSTT::EntityTopology &topo,
                                     Mesquite::MsqError &err)
  {
    switch(topo) {
    case TSTT::TRIANGLE:
      return Mesquite::TRIANGLE;
      break;
    case TSTT::QUADRILATERAL:
      return Mesquite::QUADRILATERAL;
      break;
    case TSTT::TETRAHEDRON:
      return Mesquite::TETRAHEDRON;
      break;
    case TSTT::HEXAHEDRON:
      return Mesquite::HEXAHEDRON;
      break;
    case TSTT::PRISM:
      return Mesquite::PRISM;
      break;
    case TSTT::PYRAMID:
      return Mesquite::PYRAMID;
      break;
    default:
      err.set_msg("Topology unsufficiently defined. "
                  "Cannot convert to a Mesquite Topology");
      return Mesquite::MIXED;
    }
  }

  
#define PRINT_TSTT_ERROR(tstt_err) { \
    PRINT_INFO("!!! TSTT interface ERROR caught by Mesquite -- \n!!! "); \
    PRINT_INFO(tstt_err.getNote().c_str()); \
    PRINT_INFO(tstt_err.getTrace().c_str()); \
    PRINT_INFO("\n"); \
}

}

using Mesquite::MsqError;

Mesquite::MeshTSTT::MeshTSTT(TSTT::LocalTSTTMesh& tstt_mesh,
                       Mesquite::MsqError& err) 
  : elementType(TSTT::ALL_TYPES),
    cachedVertex(NULL)
  
{
  try {
    tsttMesh = tstt_mesh;
    std::string fixed_tag("fixed");
    fixedVertexTag = tsttMesh.tagGetHandle(fixed_tag);
    std::string boundary_tag("boundary");
    boundaryVertexTag = tsttMesh.tagGetHandle(boundary_tag);
    cout << "boundaryVertexTag: " << *((int*)boundaryVertexTag) << endl; //dbg
    oneEntity = ::SIDL::array<EntityHandle>::create1d(1);
    oneTagValue = ::SIDL::array<TagHandle>::create1d(1);
    oneInt = ::SIDL::array<int32_t>::create1d(1);
    oneTopo = ::SIDL::array<TSTT::EntityTopology>::create1d(1);
    threeDoubles = ::SIDL::array<double>::create1d(3);

    cachedAdjEntArray = ::SIDL::array<EntityHandle>::create1d(0);
    adjCsrPt = ::SIDL::array<int32_t>::create1d(0);
    adjCsrDat = ::SIDL::array<int32_t>::create1d(0);

    // Associates a vertex byte flag to all vertices in the mesh.
    std::string vertex_byte_tag("MsqVtxByteTag");
    unsigned char tag_value = 0;
    tsttMesh.tagCreate(vertex_byte_tag, (int)sizeof(unsigned char),
                              &tag_value, vertexByteTag);
    SIDL::array<EntityHandle> vertices = ::SIDL::array<EntityHandle>::create1d(0);
    tsttMesh.entitysetGetEntities(ENTIRE_MESH,
                                   ::TSTT::VERTEX,
                                   ::TSTT::ALL_TOPOLOGIES,
                                   vertices);
    tsttMesh.entityAddTag(vertices, vertexByteTag);
    
  }
  catch (TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

Mesquite::MeshTSTT::~MeshTSTT() 
{
  tsttMesh.tagDelete(vertexByteTag, true);
}

    
// Returns whether this mesh lies in a 2D or 3D coordinate system.
int Mesquite::MeshTSTT::get_geometric_dimension(Mesquite::MsqError &/*err*/) const
{
  int32_t d;
  try {
    d = tsttMesh.getGeometricDimension();
  }
  catch (TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
  return (int)d;
}
    
// Returns the number of vertices.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::get_total_vertex_count" 
size_t Mesquite::MeshTSTT::get_total_vertex_count(Mesquite::MsqError &/*err*/) const
{
  int32_t nv=0;
  try {
    // gets number of vertices for whole mesh.
    nv = tsttMesh.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::VERTEX);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
    throw;
  }
  return (size_t)nv;
}

//! Returns the number of regions or number of faces if there are no regions.
//! Sets data member elementType to region or face.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::get_total_element_count" 
size_t Mesquite::MeshTSTT::get_total_element_count(Mesquite::MsqError &err) const 
{
  int32_t ne=0;
  try {
    // query nb of regions (3D elements)
    ne = tsttMesh.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::REGION);
    if (ne!=0) {
      elementType = TSTT::REGION;
    }
    // If there isn't any region, the number of elements is the number of faces (2D elements)
    else {
      ne = tsttMesh.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::FACE);
      if (ne!=0) {
        elementType = TSTT::FACE;
      }
      else
        err.set_msg("No 2D or 3D elements available in mesh.");
    }
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
  return (size_t)ne;
}
    
// Fills array with handles to all vertices in the mesh.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::get_all_vertices"
void Mesquite::MeshTSTT::get_all_vertices(
  Mesquite::Mesh::VertexHandle *vert_array,
  size_t array_size, Mesquite::MsqError &err)
{
  try {
    // Checking size compatibility.
    size_t nv = get_total_vertex_count(err); MSQ_CHKERR(err);
    if (array_size > nv)
      array_size = nv; // only copies existing vertices.
    if (array_size < nv) {
      err.set_msg("Array of insufficient size. "
                  "Returning incomplete vertex list");
    }

    // Creating borrowed SIDL array with Mesquite memory space.
    int32_t lower= 0;
    int32_t upper = array_size-1;
    int32_t stride = 1;
    ::SIDL::array<void*> vert_array_b;
    vert_array_b.borrow(vert_array, 1, &lower, &upper, &stride);
    
    // Retrieves array of vertices from TSTT interface.
    tsttMesh.entitysetGetEntities(ENTIRE_MESH,
                                   ::TSTT::VERTEX,
                                   ::TSTT::ALL_TOPOLOGIES,
                                   vert_array_b);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Fills array with handles to all elements in the mesh.
void Mesquite::MeshTSTT::get_all_elements(
  Mesquite::Mesh::ElementHandle *elem_array,
  size_t array_size, Mesquite::MsqError &err)
{
  try {
    // Checking size compatibility.
    size_t ne = get_total_element_count(err); MSQ_CHKERR(err);
    if (array_size > ne)
      array_size = ne; // only copies existing vertices.
    if (array_size < ne) {
      err.set_msg("Array of insufficient size. "
                  "Returning incomplete vertex list");
    }

    // Creating borrowed SIDL array with Mesquite memory space.
    int32_t lower= 0;
    int32_t upper = array_size-1;
    int32_t stride = 1;
    ::SIDL::array<void*> elem_array_b;
    elem_array_b.borrow(elem_array, 1, &lower, &upper, &stride);
    
    // Retrieves array of vertices from TSTT interface.
    tsttMesh.entitysetGetEntities(ENTIRE_MESH,
                                         elementType, // set by get_total_element_count
                                         ::TSTT::ALL_TOPOLOGIES,
                                         elem_array_b);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Returns a pointer to an iterator that iterates over the
// set of all vertices in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If vertices are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
Mesquite::VertexIterator* Mesquite::MeshTSTT::vertex_iterator(MsqError &/*err*/)
{
  return new MeshTSTT_EntityIterator(tsttMesh, TSTT::VERTEX);
}
    
// Returns a pointer to an iterator that iterates over the
// set of all top-level elements in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If elements are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
Mesquite::ElementIterator* Mesquite::MeshTSTT::element_iterator(MsqError &/*err*/)
{
  return new MeshTSTT_EntityIterator(tsttMesh, elementType);
}

//************ Vertex Properties ********************
// Returns true or false, indicating whether the vertex
// is allowed to be repositioned.  True indicates that the vertex
// is fixed and cannot be moved.  Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
bool Mesquite::MeshTSTT::vertex_is_fixed(Mesquite::Mesh::VertexHandle vertex, MsqError &err)
{
  try {
    int32_t tag_size;
    oneEntity.set(0,vertex);
    tsttMesh.entityGetTagData(oneEntity, fixedVertexTag,
                                     oneTagValue, tag_size);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
  return (bool)(oneTagValue.get(0));
}

// Returns true or false, indicating whether the vertex
// is on the boundary.  Boundary nodes may be treated as
// a special case by some algorithms or culling methods.
// Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
bool Mesquite::MeshTSTT::vertex_is_on_boundary(
  Mesquite::Mesh::VertexHandle vertex, MsqError &err)
{
  try {    
    int32_t tag_size = sizeof(int32_t);
    oneEntity.set(0,vertex);

    //dbg
//     EntityHandle toto = oneEntity.get(0);
//     Mesquite::Vector3D coords_2;
//     vertex_get_coordinates(toto, coords_2, err); MSQ_CHKERR(err);
//     cout << "\ncoords: " << coords_2 << endl;
//     int* v = new int(0); oneTagValue.set(0,(void*)v); 
//     tsttMesh.entitySetTagData(oneEntity, boundaryVertexTag,
//                                      oneTagValue, tag_size);
    
    tsttMesh.entityGetTagData(oneEntity, boundaryVertexTag,
                                     oneTagValue, tag_size);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
  return (bool)(*(int*)(oneTagValue.get(0)));
}

// Get/set location of a vertex
void Mesquite::MeshTSTT::vertex_get_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  Vector3D &coordinates, MsqError &err)
{
  try {
    oneEntity.set(0, vertex);
    ::TSTT::StorageOrder order = ::TSTT::INTERLEAVED;
    tsttMesh.entityGetVertexCoordinates(oneEntity,
                               order, threeDoubles);
    // Turns SIDL array into a Vector3D.
    coordinates[0] = threeDoubles.get(0);
    coordinates[1] = threeDoubles.get(1);
    coordinates[2] = threeDoubles.get(2);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

void Mesquite::MeshTSTT::vertex_set_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  const Vector3D &coordinates, MsqError &err)
{
  try {
    // Turns Vector3D into SIDL array.
    threeDoubles.set(0, coordinates[0]);
    threeDoubles.set(1, coordinates[1]);
    threeDoubles.set(2, coordinates[2]);

    oneEntity.set(0, vertex);
    tsttMesh.setVertexCoordinates(oneEntity,
                                  ::TSTT::INTERLEAVED,
                                  threeDoubles);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Each vertex has a byte-sized flag that can be used to store
// flags.  This byte's value is neither set nor used by the mesh
// implementation.  It is intended to be used by Mesquite algorithms.
// Until a vertex's byte has been explicitly set, its value is 0.
void Mesquite::MeshTSTT::vertex_set_byte (
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char byte, MsqError &err)
{
  try {
    oneEntity.set(0, vertex);
    oneTagValue.set(0, &byte);
    tsttMesh.entitySetTagData(oneEntity, vertexByteTag,
                                     oneTagValue, sizeof(unsigned char));
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

void Mesquite::MeshTSTT::vertices_set_byte (
  Mesquite::Mesh::VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
  try {
    // Creating borrowed SIDL array holding vertices.
    int32_t lower= 0;
    int32_t upper = array_size-1;
    int32_t stride = 1;
    ::SIDL::array<void*> vert_array_b;
    vert_array_b.borrow(vert_array, 1, &lower, &upper, &stride);
    
    // Creating borrowed SIDL array holding tag datas.
    ::SIDL::array<void*> byte_array_b;
    byte_array_b.borrow((void**)&byte_array, 1, &lower, &upper, &stride);

    // set tag data
    tsttMesh.entitySetTagData(vert_array_b, vertexByteTag,
                                     byte_array_b, sizeof(unsigned char));
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Retrieve the byte value for the specified vertex or vertices.
// The byte value is 0 if it has not yet been set via one of the
// *_set_byte() functions.
void Mesquite::MeshTSTT::vertex_get_byte(
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char *byte, MsqError &err)
{
  try {
    oneEntity.set(0, vertex);
    oneTagValue.set(0, &byte);
    int32_t tag_size = sizeof(unsigned char);
    tsttMesh.entityGetTagData(oneEntity, vertexByteTag,
                                     oneTagValue, tag_size);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

void Mesquite::MeshTSTT::vertices_get_byte(
  Mesquite::Mesh::VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
  try {
    // Creating borrowed SIDL array holding vertices.
    int32_t lower= 0;
    int32_t upper = array_size-1;
    int32_t stride = 1;
    ::SIDL::array<void*> vert_array_b;
    vert_array_b.borrow(vert_array, 1, &lower, &upper, &stride);
    
    // Creating borrowed SIDL array holding tag datas.
    ::SIDL::array<void*> byte_array_b;
    byte_array_b.borrow((void**)&byte_array, 1, &lower, &upper, &stride);

    // retrieve tag data
    int32_t tag_size = sizeof(unsigned char);
    tsttMesh.entityGetTagData(vert_array_b, vertexByteTag,
                                     byte_array_b, tag_size);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}


//**************** Vertex Topology *****************



// void Mesquite::MeshTSTT::create_vertex_to_element_data()
// {
//   if (v2E)
//     return;
  
//   v2eOffset = new size_t[vertexCount + 1];
//   v2E = new size_t[totalVertexUses];
  
//     // Initialize each use count to zero
//   memset(v2eOffset, 0, (vertexCount+1)*sizeof(size_t));
  
//   size_t elem_num, vert_num;
  
//     // Go through each element, keep track of how many times
//     // each vertex is used.
//   for (elem_num = elementCount; elem_num--; )
//   {
//     for (vert_num =
//            Mesquite::vertices_in_topology(elementArray[elem_num].mType);
//          vert_num--;
//          )
//     {
//       v2eOffset[elementArray[elem_num].vertexIndices[vert_num]]++;
//     }
//   }
  
//     // Convert the uses counts to array offsets.
//   elem_num = 0;
//   for (vert_num = 0; vert_num < vertexCount; vert_num++)
//   {
//     size_t temp = v2eOffset[vert_num];
//     v2eOffset[vert_num] = elem_num;
//     elem_num += temp;
//   }
//   v2eOffset[vertexCount] = totalVertexUses;

//     // Use newVertIndices to store how many elements
//     // have already been added to v2E for each vertex.
//     // Should already be initialized to zero
  
//     // Finally, store the v2E data
//   for (elem_num = 0; elem_num < elementCount; elem_num++)
//   {
//     for (vert_num =
//            Mesquite::vertices_in_topology(elementArray[elem_num].mType);
//          vert_num--;
//          )
//     {
//       size_t vert_index = elementArray[elem_num].vertexIndices[vert_num];
//       v2E[v2eOffset[vert_index] + newVertIndices[vert_index]++] = elem_num;
//     }
//   }
  
//     // Reset newVertIndices to zero
//   memset(newVertIndices, 0, vertexCount*sizeof(size_t));
// }



// Gets the number of elements attached to this vertex.
// Useful to determine how large the "elem_array" parameter
// of the vertex_get_attached_elements() function must be.
size_t Mesquite::MeshTSTT::vertex_get_attached_element_count(
  Mesquite::Mesh::VertexHandle vertex, MsqError &err) const
{
  try {
    cachedAdjEntArray = ::SIDL::array<EntityHandle>::create1d(0);
    adjCsrPt = ::SIDL::array<int32_t>::create1d(0);
    adjCsrDat = ::SIDL::array<int32_t>::create1d(0);
    oneEntity.set(0, vertex);
    tsttMesh.entityGetAdjacencies(oneEntity, elementType,
                                         cachedAdjEntArray,
                                         adjCsrPt, adjCsrDat);
    cachedVertex = vertex;
    assert(cachedAdjEntArray.dimen()==1);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }

  return static_cast<size_t>(cachedAdjEntArray.upper(0)-cachedAdjEntArray.lower(0)+1);
}

// Gets the elements attached to this vertex.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertex_get_attached_elements" 
void Mesquite::MeshTSTT::vertex_get_attached_elements(
  Mesquite::Mesh::VertexHandle vertex,
  Mesquite::Mesh::ElementHandle* elem_array,
  size_t sizeof_elem_array, Mesquite::MsqError &err)
{
  // check that adjacencies are cached
  if (vertex != cachedVertex) {
    err.set_msg("argument vertex different from latest call to "
                "vertex_get_attached_element_count. "
                "Cannot use cached data.");
    return;
  }
  // Checks that cached array will fit.
  if (sizeof_elem_array != static_cast<size_t>(cachedAdjEntArray.upper(0) - cachedAdjEntArray.lower(0)+1)) {
    err.set_msg("elem_array is not the right size. "
                "Use vertex_get_attached_element_count first.");
    return;
  }

  // fills up array from previously cached array. 
  for (size_t i=0; i<sizeof_elem_array; ++i) {
    elem_array[i] = cachedAdjEntArray.get(i);
  }
}


// Gets the number of vertices in this element.
// This data can also be found by querying the
// element's topology and getting the number
// of vertices per element for that topology type.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::element_get_attached_vertex_count" 
size_t Mesquite::MeshTSTT::element_get_attached_vertex_count(
  Mesquite::Mesh::ElementHandle elem,
  Mesquite::MsqError &err) const
{
  Mesquite::EntityTopology topo_msq;
  topo_msq = this->element_get_topology(elem, err); MSQ_CHKERR(err);
  return Mesquite::MsqMeshEntity::vertex_count(topo_msq);
}

// Returns the vertices that are part of the topological definition of each
// element in the "elem_handles" array.  When this function is called, the
// following must be true:
//   a) "elem_handles" points at an array of "num_elems" element handles.
//   b) "vert_handles" points at an array of size "sizeof_vert_handles"
//   c) "csr_data" points at an array of size "sizeof_csr_data"
//   d) "csr_offsets" points at an array of size "num_elems+1"
//      
// When this function returns, adjacency information will be stored
// in csr format:
//    a) "vert_handles" stores handles to all vertices found in one
//       or more of the elements.  Each vertex appears only
//       once in "vert_handles", even if it is in multiple elements.
//    b) "sizeof_vert_handles" is set to the number of vertex
//       handles placed into "vert_handles".
//    c) "sizeof_csr_data" is set to the total number of vertex uses (for
//       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
//       the two triangles share some vertices).
//    c) "csr_offsets" is filled such that csr_offset[i] indicates the location
//       of entity i's first adjacency in "csr_data".  The number of vertices
//       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
//       reason, csr_offsets[num_elems] is set to the new value of
//       "sizeof_csr_data".
//    d) "csr_data" stores integer offsets which give the location of
//       each adjacency in the "vert_handles" array.
//
// As an example of how to use this data, you can get the handle of the first
// vertex in element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ]
//
// and the second vertex of element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ]
// 
void Mesquite::MeshTSTT::elements_get_attached_vertices(
  Mesquite::Mesh::ElementHandle *elem_handles,
  size_t num_elems,
  Mesquite::Mesh::VertexHandle *vert_handles,
  size_t &sizeof_vert_handles,
  size_t *csr_data,
  size_t &sizeof_csr_data,
  size_t *csr_offsets,
  Mesquite::MsqError &err)
{
  if (num_elems == 0)
    return;

  try {
    // Creating borrowed SIDL arrays with arrays passed as arguments.
    int32_t lower= 0;
    int32_t stride = 1;
    int32_t upper = num_elems-1;
    ::SIDL::array<void*> elem_handles_b;
    elem_handles_b.borrow(elem_handles, 1, &lower, &upper, &stride);
    upper = sizeof_vert_handles-1;
    ::SIDL::array<void*> vert_handles_s;
    vert_handles_s = ::SIDL::array<void*>::create1d(0);
    upper = num_elems;
    ::SIDL::array<int32_t> csr_offsets_b;
    csr_offsets_b.borrow((int32_t*)csr_offsets, 1, &lower, &upper, &stride);
    ::SIDL::array<int32_t> dummy;
    dummy = ::SIDL::array<int32_t>::create1d(0);

    tsttMesh.entityGetAdjacencies(elem_handles_b, ::TSTT::VERTEX,
                                         vert_handles_s, csr_offsets_b,
                                         dummy);

    // TODO : assert csr_offsets_b has not been reallocated by TSTT implementation.

    // TODO: Below, distance() is costly. There are faster ways ... 
     
    // Converts the TSTT 2 arrays format with repeated adhjacency handles into the
    // mesquite 3 arrays CSR format without repeated handles.
    size_t d=0;
    std::set<Mesquite::Mesh::VertexHandle> vert_handles_csr;
    for (size_t e=0; e<num_elems; ++e) {
      for (size_t v=csr_offsets[e]; v<csr_offsets[e+1]; ++v) { 

	// checking we're not writing beyond the index array 
	if(d>sizeof_csr_data) {
	  err.set_msg("insuficient size for csr_data array");
	  return;
	}
	
	// std::set::insert() returns a pair with the location of the entry
	// and whether an insert was actually done.
	std::pair<std::set<Mesquite::Mesh::VertexHandle>::iterator,bool> status;
        status = vert_handles_csr.insert(vert_handles_s[v]);
        csr_data[d]=distance(vert_handles_csr.begin(),
                             status.first);

	// if a new unique handle has been inserted
	// increment the indexes equal or bigger in the index array
	if (status.second == true)
          for (size_t j=0; j<d; ++j) 
	    if (csr_data[j]>=csr_data[d]) 
	      csr_data[j]+=1;
			      
        ++d;
      }
    }
    sizeof_csr_data=d;
     
    // Now just copies exactly the std::set into the array 
    // given as an (output) argument.
    size_t v=0;
    std::set<Mesquite::Mesh::VertexHandle>::iterator vertices;
    for (vertices=vert_handles_csr.begin();
        vertices!=vert_handles_csr.end();
         ++vertices) {
      vert_handles[v] = *vertices;
      assert(v<=sizeof_vert_handles);
      ++v;
    }
    assert(v!=0);
    sizeof_vert_handles=v;

  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Identifies the vertices attached to this element by returning
// each vertex's global index.  The vertex's global index indicates
// where that element can be found in the array returned by
// Mesh::get_all_vertices.
void Mesquite::MeshTSTT::element_get_attached_vertex_indices(
  Mesquite::Mesh::ElementHandle element,
  size_t *index_array,
  size_t array_size, MsqError &err)
{
  try {
    int32_t lower= 0;
    int32_t stride = 1;
    int32_t upper = Mesquite::MSQ_MAX_NUM_VERT_PER_ENT -1;
    ::SIDL::array<int32_t> index_array_b;
    index_array_b.borrow((int32_t*)index_array, 1, &lower, &upper, &stride);
    oneEntity.set(0, element);

    //! \todo : THIS IS A HACK : the first argument is not valid for a real TSTT interface.
    //! An EntitySet with a single member element should be used instead.
    //! This hack catters to a temporary AOMD implementation of TSTT without EntitySet
    tsttMesh.entitysetGetEntityVertexCoordinateIndices(
                                element, elementType, ::TSTT::ALL_TOPOLOGIES,
                                oneInt, index_array_b,
                                oneTopo);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Returns the topology of the given entity.
Mesquite::EntityTopology Mesquite::MeshTSTT::element_get_topology(
  Mesquite::Mesh::ElementHandle element, Mesquite::MsqError &err) const
{
  Mesquite::EntityTopology topo_msq;
  try {
    oneEntity.set(0, element);
    tsttMesh.entityGetTopology(oneEntity, oneInt);
    TSTT::EntityTopology topo_tstt = (TSTT::EntityTopology)oneInt.get(0);
    topo_msq = mesquite_equivalent_topology(topo_tstt, err); MSQ_CHKERR(err);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
  return topo_msq;
}

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
void Mesquite::MeshTSTT::elements_get_topologies(
  Mesquite::Mesh::ElementHandle *element_handle_array,
  Mesquite::EntityTopology *element_topologies,
  size_t num_elements, MsqError &err)
{
//   for (size_t i = 0; i < num_elements; i++)
//   {
//     element_topologies[i] =
//       reinterpret_cast<MeshTSTT::Element*>(element_handle_array[i])->mType;
//   }
}

//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void Mesquite::MeshTSTT::release_entity_handles(
  Mesquite::Mesh::EntityHandle */*handle_array*/,
  size_t /*num_handles*/, MsqError &err)
{
    // Do nothing
}

// Instead of deleting a Mesh when you think you are done,
// call release().  In simple cases, the implementation could
// just call the destructor.  More sophisticated implementations
// may want to keep the Mesh object to live longer than Mesquite
// is using it.
void Mesquite::MeshTSTT::release()
{
  tsttMesh.deleteRef();
//  delete this;
}
