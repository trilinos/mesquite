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


  \author Thomas Leurent
  \date   2003-06-12
*/

//kkc 040203 - changed implementation to use more general TSTT interfaces

#include <set>
//#include <pair>

#include "MeshTSTT.hpp"
//kkc 040203 now included in MeshTSTT.hpp #include "TSTT.hh"

#include "MsqTimer.hpp"

using std::cout;
using std::cerr;
using std::endl;

namespace // unnamed namespace (scope is the file only)
{

#define ENTIRE_MESH 0  // passing the null pointer as the first argument of
                   // TSTT::entitySet functions queries the data for the whole mesh

#define PRINT_TSTT_ERROR(tstt_err) { \
    Message::print_info("!!! TSTT interface ERROR caught by Mesquite -- \n!!! "); \
    Message::print_info("%s\n",tstt_err.getNote().c_str()); \
    std::string errs; tstt_err.getErrorDescription(errs); \
    Message::print_info("%s\n",errs.c_str()); \
    Message::print_info("%s\n",tstt_err.getTrace().c_str()); \
}

  
  class MeshTSTT_EntityIterator : public Mesquite::EntityIterator
  {
  public:
    MeshTSTT_EntityIterator(TSTT::Mesh& tstt_mesh,
                            const TSTT::EntityType& entity_type)
      : iteratorEntityType(entity_type)
    {
      try {
        iteratorMesh = tstt_mesh; 
	if ( !iteratorMesh )
	  {
	    //kkc XXX is this how errors should be handled?
	    TSTT::Error err = TSTT::Error::_create();
	    err.set(1,"MeshTSTT_EntityIterator::MeshTSTT_EntityIterator : TSTT::Mesh does not"
		    " implement TSTT::AdvancedEntitySetQuery");
	    throw err;
	  }

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
      catch (TSTT::Error &tstt_err) {
        PRINT_TSTT_ERROR(tstt_err);
      }
    }
    
      // Moves the iterator back to the first
      // entity in the list.
    virtual void restart()
      {
        try{
          iteratorMesh.entitysetDestroyWorksetIterator(
                                 tsttWorksetIterator);
          iteratorMesh.entitysetInitializeWorksetIterator(
                                         ENTIRE_MESH,
                                         iteratorEntityType,
                                         ::TSTT::ALL_TOPOLOGIES,
                                         1, // requested workset size
                                         tsttWorksetIterator);
        }
        catch (TSTT::Error &tstt_err) {
          PRINT_TSTT_ERROR(tstt_err);
        }
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
        try {
          bool more = iteratorMesh.entitysetGetNextWorkset(tsttWorksetIterator, entityHandle);
          if (!more)
            entityHandle.set(0, NULL);
        }
        catch (TSTT::Error &tstt_err) {
          PRINT_TSTT_ERROR(tstt_err);
        }
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
    //kkc 040203    ::TSTT::LocalTSTTMesh iteratorMesh;
    ::TSTT::AdvancedEntitySetQuery iteratorMesh;
    ::TSTT::EntityType iteratorEntityType;
    mutable ::SIDL::array<void*> entityHandle; // Array has one entry only. That entry is
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

} // end of unnamed namespace (scope is the file only)

using Mesquite::MsqError;

#undef __FUNC__
#define __FUNC__ "MeshTSTT::MeshTSTT" 
Mesquite::MeshTSTT::MeshTSTT(TSTT::Mesh& tstt_mesh,
                       Mesquite::MsqError& err) 
  : elementType(TSTT::ALL_TYPES),
    cachedVertex(NULL)
  
{
  try {
    //kkc 040203 cast the appropriate mesh interfaces, make sure they are implemented
    //           by checking the result of the cast
    //kkc tsttMesh = tstt_mesh;

    this->tsttMesh = tstt_mesh;
    tsttCoreQuery = tsttMesh;
    if ( !tsttCoreQuery ) 
      {
	TSTT::Error err = TSTT::Error::_create();
	err.set(1,"Mesquite::MeshTSTT::MeshTSTT : the TSTT::Mesh does not implement TSTT::CoreEntitySetQuery");
	throw err;
      }

    tsttAdvQuery  = tsttMesh;
    if ( !tsttAdvQuery ) 
      {
	TSTT::Error err = TSTT::Error::_create();
	err.set(1,"Mesquite::MeshTSTT::MeshTSTT : the TSTT::Mesh does not implement TSTT::AdvancedEntitySetQuery");
	throw err;
      }

    tsttTag       = tsttMesh;
    if ( !tsttTag ) 
      {
	TSTT::Error err = TSTT::Error::_create();
	err.set(1,"Mesquite::MeshTSTT::MeshTSTT : the TSTT::Mesh does not implement TSTT::Tag");
	throw err;
      }

    tsttModMesh   = tsttMesh;
    if ( !tsttModMesh ) 
      {
	TSTT::Error err = TSTT::Error::_create();
	err.set(1,"Mesquite::MeshTSTT::MeshTSTT : the TSTT::Mesh does not implement TSTT::ModifiableMesh");
	throw err;
      }

    //kkc 040203 if we get here then presumably all the interfaces we need are implemented by the
    //           TSTT::Mesh implementation we have been given.  

    std::string fixed_tag("fixed");
    //kkc 040203    fixedVertexTag = tsttMesh.tagGetHandle(fixed_tag);
    //kkc 040203 not used?    fixedVertexTag = tsttTag.tagGetHandle(fixed_tag);

    std::string boundary_tag("boundary");
    //kkc 040203    boundaryVertexTag = tsttMesh.tagGetHandle(boundary_tag);
    boundaryVertexTag = tsttTag.tagGetHandle(boundary_tag);

//    cout << "boundaryVertexTag: " << *((int*)boundaryVertexTag) << endl; //dbg
    oneEntity = ::SIDL::array<EntityHandle>::create1d(1);
    oneTagValue = ::SIDL::array<TagHandle>::create1d(1);
    oneInt = ::SIDL::array<int32_t>::create1d(1);
    oneTopo = ::SIDL::array<TSTT::EntityTopology>::create1d(1);
    threeDoubles = ::SIDL::array<double>::create1d(3);

    //kkc actually set the following to uninitialized arrays
    //cachedAdjEntArray = ::SIDL::array<EntityHandle>::create1d(0);
    //adjCsrPt = ::SIDL::array<int32_t>::create1d(0);
    //adjCsrDat = ::SIDL::array<int32_t>::create1d(0);
    cachedAdjEntArray = ::SIDL::array<EntityHandle>();
    adjCsrPt = ::SIDL::array<int32_t>();
    adjCsrDat = ::SIDL::array<int32_t>();

    // Associates a vertex byte flag to all vertices in the mesh.
    std::string vertex_byte_tag("MsqVtxByteTag");
    unsigned char tag_value = 0;
    //kkc 040203    tsttMesh.tagCreate(vertex_byte_tag, (int)sizeof(unsigned char),
    //                              &tag_value, vertexByteTag);
    tsttTag.tagCreate(vertex_byte_tag, (int)sizeof(unsigned char),
                              &tag_value, vertexByteTag);
    //    SIDL::array<EntityHandle> vertices = ::SIDL::array<EntityHandle>::create1d(0);
    SIDL::array<EntityHandle> vertices;//kkc = ::SIDL::array<EntityHandle>::create1d(0);

    //kkc 040203    tsttMesh.entitysetGetEntities(ENTIRE_MESH,
    //                                   ::TSTT::VERTEX,
    //                                   ::TSTT::ALL_TOPOLOGIES,
    //                                   vertices);
    tsttCoreQuery.entitysetGetEntities(ENTIRE_MESH,
					::TSTT::VERTEX,
					::TSTT::ALL_TOPOLOGIES,
					vertices);

    //kkc 040203    tsttMesh.entityAddTag(vertices, vertexByteTag);
    tsttTag.entityAddTag(vertices, vertexByteTag);

    int nb_vertices = vertices.upper(0) + 1;
    mZeros = new unsigned char[nb_vertices];
    for (int i=0; i<nb_vertices; ++i)
      mZeros[i] = 0;
    ::SIDL::array<void*> vtx_byte_values = ::SIDL::array<void*>::create1d(nb_vertices);
    for (int i=0; i<nb_vertices; ++i) {
      vtx_byte_values.set(i, &mZeros[i]);
    }
    //kkc 040203    tsttMesh.entitySetTagData(vertices, vertexByteTag, vtx_byte_values,
    //                              (int32_t)sizeof(unsigned char));
    tsttTag.entitySetTagData(vertices, vertexByteTag, vtx_byte_values,
                              (int32_t)sizeof(unsigned char));

    nbElements = get_total_element_count(err); MSQ_CHKERR(err);
  }
  catch (TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

#undef __FUNC__
#define __FUNC__ "MeshTSTT::~MeshTSTT" 
Mesquite::MeshTSTT::~MeshTSTT() 
{
  //kkc 040203  tsttMesh.tagDelete(vertexByteTag, true);
  tsttTag.tagDelete(vertexByteTag, true);
}

    
// Returns whether this mesh lies in a 2D or 3D coordinate system.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::get_geometric_dimension" 
int Mesquite::MeshTSTT::get_geometric_dimension(Mesquite::MsqError &/*err*/) const
{
  int32_t d=0;
  try {
    //kkc 040203    d = tsttMesh.getGeometricDimension();
    d = tsttCoreQuery.getGeometricDimension();
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
    //kkc 040203    nv = tsttMesh.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::VERTEX);
    nv = tsttAdvQuery.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::VERTEX);
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
    //kkc 040203    ne = tsttMesh.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::REGION);
    ne = tsttAdvQuery.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::REGION);
    if (ne!=0) {
      elementType = TSTT::REGION;
    }
    // If there isn't any region, the number of elements is the number of faces (2D elements)
    else {
      //kkc 040203      ne = tsttMesh.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::FACE);
      ne = tsttAdvQuery.entitysetGetNumberEntityOfType(ENTIRE_MESH, ::TSTT::FACE);
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
    //kkc 040203    tsttMesh.entitysetGetEntities(ENTIRE_MESH,
    //                                   ::TSTT::VERTEX,
    //                                   ::TSTT::ALL_TOPOLOGIES,
    //                                   vert_array_b);
    tsttCoreQuery.entitysetGetEntities(ENTIRE_MESH,
					::TSTT::VERTEX,
					::TSTT::ALL_TOPOLOGIES,
					vert_array_b);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Fills array with handles to all elements in the mesh.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::get_all_all_elements"
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
    //kkc 040203    tsttMesh.entitysetGetEntities(ENTIRE_MESH,
    //                                         elementType, // set by get_total_element_count
    //                                         ::TSTT::ALL_TOPOLOGIES,
    //                                         elem_array_b);
    tsttCoreQuery.entitysetGetEntities(ENTIRE_MESH,
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
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertex_iterator"
Mesquite::VertexIterator* Mesquite::MeshTSTT::vertex_iterator(MsqError &/*err*/)
{
  //kkc 040203  return new MeshTSTT_EntityIterator(tsttMesh, TSTT::VERTEX);
  return new MeshTSTT_EntityIterator(tsttMesh, TSTT::VERTEX);
}
    
// Returns a pointer to an iterator that iterates over the
// set of all top-level elements in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If elements are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::element_iterator"
Mesquite::ElementIterator* Mesquite::MeshTSTT::element_iterator(MsqError &/*err*/)
{
  //kkc 040203  return new MeshTSTT_EntityIterator(tsttMesh, elementType);
  return new MeshTSTT_EntityIterator(tsttMesh, elementType);
}

//************ Vertex Properties ********************
// Returns true or false, indicating whether the vertex
// is allowed to be repositioned.  True indicates that the vertex
// is fixed and cannot be moved.  Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertex_is_fixed"
bool Mesquite::MeshTSTT::vertex_is_fixed(Mesquite::Mesh::VertexHandle vertex, MsqError &err)
{
//   try {
//     int32_t tag_size;
//     oneEntity.set(0,vertex);
//     tsttMesh.entityGetTagData(oneEntity, fixedVertexTag,
//                                      oneTagValue, tag_size);
//   }
//   catch(::TSTT::Error &tstt_err) {
//     PRINT_TSTT_ERROR(tstt_err);
//   }
//   return (bool)(oneTagValue.get(0));

  return false; //! \todo vertex_is_fixed is not implemented for now. Not really used either.
}

// Returns true or false, indicating whether the vertex
// is on the boundary.  Boundary nodes may be treated as
// a special case by some algorithms or culling methods.
// Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertices_are_on_boundary"
void Mesquite::MeshTSTT::vertices_are_on_boundary(
  Mesquite::Mesh::VertexHandle vert_array[], bool on_bnd[],
  size_t num_vtx, MsqError &err)
{
  try {
    FUNCTION_TIMER_START(__FUNC__);

    int32_t lower=0, upper=num_vtx-1, stride=1;
    ::SIDL::array<void*> vert_array_b;
    vert_array_b.borrow(vert_array, 1, &lower, &upper, &stride);
    
    ::SIDL::array<void*> tag_s = ::SIDL::array<void*>::create1d(num_vtx);
    int32_t tag_size = sizeof(int32_t);

    //dbg
//     EntityHandle toto = oneEntity.get(0);
//     Mesquite::Vector3D coords_2;
//     vertex_get_coordinates(toto, coords_2, err); MSQ_CHKERR(err);
//     cout << "\ncoords: " << coords_2 << endl;
//     int* v = new int(0); oneTagValue.set(0,(void*)v); 
//    tsttMesh.entitySetTagData(oneEntity, boundaryVertexTag,
//                                     oneTagValue, tag_size);
    
    tsttTag.entityGetTagData(vert_array_b, boundaryVertexTag,
                              tag_s, tag_size);

    for (size_t i=0; i<num_vtx; ++i)
      on_bnd[i] = (bool)(*(int*)(tag_s.get(i)));

    FUNCTION_TIMER_END();
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Get/set location of a vertex
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertices_get_coordinates"
void Mesquite::MeshTSTT::vertices_get_coordinates(
  Mesquite::Mesh::VertexHandle vert_array[],
  MsqVertex* const &coordinates, const size_t &num_vtx, MsqError &err)
{
  try {
    FUNCTION_TIMER_START(__FUNC__);
    int32_t lower=0, upper=num_vtx-1, stride=1;
    ::SIDL::array<void*> vert_array_b;
    vert_array_b.borrow(vert_array, 1, &lower, &upper, &stride);

    ::SIDL::array<double> coords_s = ::SIDL::array<double>::create1d(num_vtx*3);
    
    ::TSTT::StorageOrder order = ::TSTT::INTERLEAVED;
    tsttCoreQuery.entityGetVertexCoordinates(vert_array_b,
                                        order, coords_s);
    
    // Turns SIDL array into a Vector3D.
    int geom_dim = get_geometric_dimension(err); MSQ_CHKERR(err);
//    cout << "Geometric Dimentsion : " << geom_dim << endl;
    for (size_t i=0; i<num_vtx; ++i) {
      coordinates[i][0] = coords_s.get(geom_dim*i);
      coordinates[i][1] = coords_s.get(geom_dim*i+1);
      if ( geom_dim==3 )
	coordinates[i][2] = coords_s.get(geom_dim*i+2);
      else
	coordinates[i][2] = 0;
    }
    FUNCTION_TIMER_END();

  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertex_set_coordinates"
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
    tsttModMesh.setVertexCoordinates(oneEntity,
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
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertex_set_byte"
void Mesquite::MeshTSTT::vertex_set_byte (
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char byte, MsqError &err)
{
  try {
    oneEntity.set(0, vertex);
    unsigned char* byte_s = new unsigned char; // leaky . may catter to AOMD implementation only
    *byte_s = byte;
    oneTagValue.set(0, byte_s);
    tsttTag.entitySetTagData(oneEntity, vertexByteTag,
			      oneTagValue, sizeof(unsigned char));
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertices_set_byte"
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
    ::SIDL::array<void*> byte_array_s = ::SIDL::array<void*>::create1d(array_size);
    for (size_t i=0; i<array_size; ++i)
      byte_array_s.set(i, &(byte_array[i]));
    
    // set tag data
    //kkc 040203    tsttMesh.entitySetTagData(vert_array_b, vertexByteTag,
    //                                     byte_array_s, sizeof(unsigned char));
    tsttTag.entitySetTagData(vert_array_b, vertexByteTag,
			      byte_array_s, sizeof(unsigned char));
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Retrieve the byte value for the specified vertex or vertices.
// The byte value is 0 if it has not yet been set via one of the
// *_set_byte() functions.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertex_get_byte"
void Mesquite::MeshTSTT::vertex_get_byte(
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char *byte, MsqError &err)
{
  try {
    oneEntity.set(0, vertex);
    int32_t tag_size = sizeof(unsigned char);
    //kkc 040203    tsttMesh.entityGetTagData(oneEntity, vertexByteTag,
    //                                     oneTagValue, tag_size);
    tsttTag.entityGetTagData(oneEntity, vertexByteTag,
			      oneTagValue, tag_size);
    *byte = *(unsigned char*)oneTagValue.get(0);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertices_get_byte"
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
    
    // Creating SIDL array holding tag datas.
    ::SIDL::array<void*> byte_array_s = ::SIDL::array<void*>::create1d(array_size);

    // retrieve tag data
    int32_t tag_size = sizeof(unsigned char);
    //kkc 040203    tsttMesh.entityGetTagData(vert_array_b, vertexByteTag,
    //                                     byte_array_s, tag_size);
    tsttTag.entityGetTagData(vert_array_b, vertexByteTag,
                                     byte_array_s, tag_size);

    for (size_t i=0; i<array_size; ++i)
      byte_array[i] = *(unsigned char*)byte_array_s.get(i);
      
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}


//**************** Vertex Topology *****************


// Gets the number of elements attached to this vertex.
// Useful to determine how large the "elem_array" parameter
// of the vertex_get_attached_elements() function must be.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::vertex_get_attached_element_count"
size_t Mesquite::MeshTSTT::vertex_get_attached_element_count(
  Mesquite::Mesh::VertexHandle vertex, MsqError &err) const
{
  try {
    //kkc 040204 set these arrays to uninitialized values
    //    cachedAdjEntArray = ::SIDL::array<EntityHandle>::create1d(0);
    //adjCsrPt = ::SIDL::array<int32_t>::create1d(0);
    //adjCsrDat = ::SIDL::array<int32_t>::create1d(0);
    cachedAdjEntArray = ::SIDL::array<EntityHandle>();
    adjCsrPt = ::SIDL::array<int32_t>();
    adjCsrDat = ::SIDL::array<int32_t>();

    oneEntity.set(0, vertex);
    //kkc 040203    tsttMesh.entityGetAdjacencies(oneEntity, elementType,
    //                                         cachedAdjEntArray,
    //                                         adjCsrPt, adjCsrDat);
    tsttAdvQuery.entityGetAdjacencies(oneEntity, elementType,
				       cachedAdjEntArray,
				       adjCsrPt);//, adjCsrDat);
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
#undef __FUNC__
#define __FUNC__ "MeshTSTT::elements_get_attached_vertices"
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
    FUNCTION_TIMER_START(__FUNC__);
    // Creating borrowed SIDL arrays with arrays passed as arguments.
    int32_t lower= 0;
    int32_t stride = 1;
    int32_t upper = num_elems-1;
    ::SIDL::array<void*> elem_handles_b;
    elem_handles_b.borrow(elem_handles, 1, &lower, &upper, &stride);
    upper = sizeof_vert_handles-1;
    ::SIDL::array<void*> vert_handles_s;
    //kkc 040204 set the array to be uninitialized vert_handles_s = ::SIDL::array<void*>::create1d(0);
    vert_handles_s = ::SIDL::array<void*>();

    upper = num_elems;
    ::SIDL::array<int32_t> csr_offsets_b;
    csr_offsets_b.borrow((int32_t*)csr_offsets, 1, &lower, &upper, &stride);
    ::SIDL::array<int32_t> dummy;
    //kkc 040204 not needed    dummy = ::SIDL::array<int32_t>::create1d(0);

    //kkc 040203    tsttMesh.entityGetAdjacencies(elem_handles_b, ::TSTT::VERTEX,
    //                                         vert_handles_s, csr_offsets_b,
    //                                         dummy);
    tsttAdvQuery.entityGetAdjacencies(elem_handles_b, ::TSTT::VERTEX,
				       vert_handles_s, csr_offsets_b);//,
    //				       dummy);

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
     FUNCTION_TIMER_END();
	  return;
	}
	
	// std::set::insert() returns a pair with the location of the entry
	// and whether an insert was actually done.
	std::pair<std::set<Mesquite::Mesh::VertexHandle>::iterator,bool> status;
	//kkc 040203        status = vert_handles_csr.insert(vert_handles_s[v]);
        status = vert_handles_csr.insert(vert_handles_s.get(v));

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
    FUNCTION_TIMER_END();
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Identifies the vertices attached to the elements by returning
// each vertex's global index.  The vertex's global index indicates
// where that element can be found in the array returned by
// Mesh::get_all_vertices. The indexes can be repeated.
// \param elems Array of element handles.
// \param num_elems number of elements in the array elems
// \param index_array Array containing the indexes of the elements vertices.
//                    Indexes can be repeated.
// \param index_array_size indicates the size of index_array
// \param offsets An array of offsets into the index_array that indicates where
//                the indexes corresponding to an element start. First entry is 0.
// For example element, to get the vertex indices of elems[3], we look at the entries
// index_array[offsets[3]] to index_array[offsets[4]] . 
#undef __FUNC__
#define __FUNC__ "MeshTSTT::elements_get_attached_vertex_indices"
void Mesquite::MeshTSTT::elements_get_attached_vertex_indices(
  Mesquite::Mesh::ElementHandle elems[],
  size_t num_elems,
  size_t index_array[],
  size_t array_size,
  size_t* offsets,
  MsqError &err)
{
  try {
    FUNCTION_TIMER_START(__FUNC__);
    int32_t zero= 0;
    int32_t upper = array_size-1;
    int32_t stride = 1;
    ::SIDL::array<int32_t> index_array_b;
    index_array_b.borrow((int32_t*)index_array, 1, &zero, &upper, &stride);

    int32_t num_elems_m1 = num_elems - 1;
    ::SIDL::array<void*> elems_b;
    elems_b.borrow(elems, 1, &zero, &num_elems_m1, &stride);
    upper = num_elems;
    ::SIDL::array<int32_t> offsets_b;
    offsets_b.borrow((int32_t*)offsets, 1, &zero, &upper , &stride);
    ::SIDL::array<TSTT::EntityTopology> topos_s = ::SIDL::array<TSTT::EntityTopology>::create1d(0);

    //! \todo Here one should create an entityset with the elements handles. !!!!!!
    //!       and use it instead of ENTIRE_MESH in the following functions call. !!!!

    assert( num_elems==nbElements );

    tsttCoreQuery.entitysetGetEntityVertexCoordinateIndices(
                          ENTIRE_MESH, elementType, ::TSTT::ALL_TOPOLOGIES,
                          offsets_b, index_array_b,
                          topos_s);
    //dbg
//     cout << "array_size: "<<array_size << endl;
//     for (int i=0; i<array_size; ++i) 
//     {
//        cout << "index_array_b: " << index_array_b.get(i) << endl;
//        cout << "index_array: " << index_array[i] << endl;	  
//     }
     
    FUNCTION_TIMER_END();
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

// Returns the topology of the given entity.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::element_get_topology"
Mesquite::EntityTopology Mesquite::MeshTSTT::element_get_topology(
  Mesquite::Mesh::ElementHandle element, Mesquite::MsqError &err) const
{
  Mesquite::EntityTopology topo_msq=MIXED;
  try {
    oneEntity.set(0, element);
    //kkc 040203    tsttMesh.entityGetTopology(oneEntity, oneInt);
    tsttAdvQuery.entityGetTopology(oneEntity, oneTopo);
    //kkc 040203    TSTT::EntityTopology topo_tstt = (TSTT::EntityTopology)oneInt.get(0);
    TSTT::EntityTopology topo_tstt = (TSTT::EntityTopology)oneTopo.get(0);
    topo_msq = mesquite_equivalent_topology(topo_tstt, err); MSQ_CHKERR(err);
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
  return topo_msq;
}

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
#undef __FUNC__
#define __FUNC__ "MeshTSTT::elements_get_topologies"
void Mesquite::MeshTSTT::elements_get_topologies(
  Mesquite::Mesh::ElementHandle *element_handle_array,
  Mesquite::EntityTopology *element_topologies,
  size_t num_elements, MsqError &err)
{
  try {
    FUNCTION_TIMER_START(__FUNC__);
    ::SIDL::array< ::TSTT::EntityTopology > topos_s;
    topos_s = ::SIDL::array< ::TSTT::EntityTopology >::create1d(num_elements);

    int32_t stride = 1;
    int32_t lower = 0;
    int32_t upper = num_elements-1;
    ::SIDL::array<void*> element_handles_b;
    element_handles_b.borrow((void**)element_handle_array, 1, &lower, &upper, &stride);

    tsttAdvQuery.entityGetTopology( element_handles_b, topos_s);
    for(size_t i=0; i<num_elements; ++i) {
      TSTT::EntityTopology topo_tstt = (TSTT::EntityTopology)topos_s.get(i);
      element_topologies[i]= mesquite_equivalent_topology(topo_tstt, err);
      MSQ_CHKERR(err);
    }
    FUNCTION_TIMER_END();
  }
  catch(::TSTT::Error &tstt_err) {
    PRINT_TSTT_ERROR(tstt_err);
  }
}

//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
#undef __FUNC__
#define __FUNC__ "MeshTSTT::release_entity_handles"
void Mesquite::MeshTSTT::release_entity_handles(
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
void Mesquite::MeshTSTT::release()
{
  // remove byte tag from all vertices. 
  //kkc 040204 the initialization is not needed
  //         SIDL::array<EntityHandle> vertices = ::SIDL::array<EntityHandle>::create1d(0);
  SIDL::array<EntityHandle> vertices;//kkc = ::SIDL::array<EntityHandle>::create1d(0);

  //kkc 040203  tsttMesh.entitysetGetEntities(ENTIRE_MESH,
  //                                ::TSTT::VERTEX,
  //                                ::TSTT::ALL_TOPOLOGIES,
  //                                vertices);
  //kkc 040203  tsttMesh.entityRemoveTag(vertices,vertexByteTag);
  tsttCoreQuery.entitysetGetEntities(ENTIRE_MESH,
                                ::TSTT::VERTEX,
                                ::TSTT::ALL_TOPOLOGIES,
                                vertices);
  tsttTag.entityRemoveTag(vertices,vertexByteTag);
//  tsttMesh.tagDelete(vertexByteTag,true); // should do the same as all of the above

//kkc 040203  tsttMesh.deleteRef();
  tsttModMesh.deleteRef();  
  tsttTag.deleteRef();  
  tsttAdvQuery.deleteRef();  
  tsttCoreQuery.deleteRef();
  tsttMesh.deleteRef();  

//  delete this;
}
