#ifndef TSTT_INTERFACE_BASE_H
#define TSTT_INTERFACE_BASE_H

namespace TSTT
{
//====================================================================
// Enumerator EntityType
//====================================================================

/**
TSTT supports zero-, one-, two-, and three-dimensional entities associated 
with a mesh infrastructure, we allow users to access these dimensional 
entities using the enumerated type EntityType which contains
  @param VERTEX a zero-dimensional entity (value=0)
  @param EDGE a one-dimensional entity  (value=1)
  @param FACE a two-dimensional entity  (value=2)
  @param REGION a three-dimensional entity  (value=3)
  @param NUMBER_OF_ENTITY_TYPES the number of entities (value=4)
  */

  enum EntityType 
  {
    VERTEX=            0,
    EDGE=              1,
    FACE=              2,
    REGION=            3,
    NUMBER_OF_ENTITY_TYPES=  4
   };

/**

An enumeration of topological TSTT entities.  Note that not all TSTT meshes 
need to support all of these element types, but they do need to be able to answer 
questions such as: 
   how many elements of each type do you have?
   and return the appropriate type for entity handles

So far we agreed to enumerate the following types - more could be 
defined as time progresses.
  @param POINT a general zero-dimensional entity (value=5)
  @param LINE a general one-dimensional entity (value=6)
  @param POLYGON a general two-dimensional element (value=7)
  @param TRIANGLE a three-sided, two-dimensional element (value=8)
  @param QUADRILATERAL a four-sided, two-dimensional element (value=9)
  @param POLYHEDRAL a general three-dimensional element (value=10)
  @param TETRAHDRON a four-sided, three-dimensional element whose faces are triangles (value=11)
  @param HEXAHEDRON a six-sided, three-dimensional element whose faces are quadrilaterals (value=12)
  @param PRISM a five-sided, three-dimensional element which has three quadrilateral faces and two triangular faces (value=13)
  @param PYRAMID a five-sided, three-dimensional element which has one quadrilateral face and four triangular faces (value=14)
  @param SEPTAHEDRON a seven-sided, three-dimensional element which has one hexagonal face and six triangular faces (value=15)
  @param UNDEFINED the value that is returned when the entity topology is unknown

*/

  enum EntityTopology
  {
    POINT=             5,
    LINE=              6,
    POLYGON=           7,
    TRIANGLE=          8, 
    QUADRILATERAL=     9, 
    POLYHEDRON=        10,
    TETRAHEDRON=       11, 
    HEXAHEDRON=        12, 
    PRISM=             13, 
    PYRAMID=           14, 
    SEPTAHEDRON=       15,
    UNDEFINED=         16
  };


//====================================================================
// Enumerator StorageOrder
//====================================================================


/** An enumerated type giving the storage order for
    arrays of data (such as coordinate information).

  @param BLOCKED for storage orders such as xxxx...yyyy....zzzz(value=0)
  @param INTERLEAVED for storage orders such as xyzxyzxyzxyz... (value=1)
  @param UNDETERMINED returned if the storage order cannot be determined (value=2)

*/
  enum StorageOrder
  {
    BLOCKED=           0,
    INTERLEAVED=       1,
    UNDETERMINED=      2
  };

//====================================================================
// Typedef Int
//====================================================================

/**
A generic typedef for integers.
 */
  typedef int Int;


//====================================================================
// Typedef Float
//====================================================================

/**
A generic typedef for floats.
 */
  typedef float Float;

//====================================================================
// Typedef Double
//====================================================================

/**
A generic typedef for doubles.
 */
  typedef double Double;


//====================================================================
// Typedef Entity_Handle
//====================================================================
  
/**
An opaque object containing enough information for the underlying implementation
to determine a unique entity.
 */
  typedef void * Entity_Handle;


//====================================================================
// Typedef cEntity_Handle
//====================================================================

/**
An const opaque object containing enough information for the underlying 
implementation to determine a unique entity.  This is used to distinquish
which subroutines can change the opaque object and which cannot.
 */
  typedef const void * cEntity_Handle;


//====================================================================
// Typedef Mesh_Handle
//====================================================================

/**
An opaque object containing enough information for the underlying implementation
to determine a unique mesh.
 */
  typedef void * Mesh_Handle;


//====================================================================
// Typedef cMesh_Handle
//====================================================================

/**
An const opaque object containing enough information for the underlying 
implementation to determine a unique mesh.  This is used to distinquish
which subroutines can change the opaque object and which cannot.
 */
  typedef const void * cMesh_Handle;


//====================================================================
// Typedef Identifier
//====================================================================

/** 
An opaque object that can be used to define unique IDs for entities.
It still needs to be decided if IDs are unique among each entity type 
(ie. can we have vertex 1 and edge 1?) or if they are unique across
all entities in the mesh.
*/
 /*typedef void * Identifier;*/


//====================================================================
// Typedef cIdentifier
//====================================================================

/** 
An const opaque object that can be used to define unique IDs for entities.
*/
 /*typedef const void * cIdentifier;*/


//====================================================================
// Typedef Entity_Workset_Iterator
//====================================================================

/**
An opaque object containing enough information for the underlying
implementation to determine a unique entity.  
 */
  typedef void * Entity_Workset_Iterator;

//====================================================================
// Typedef cEntity_Workset_Iterator
//====================================================================

/**
An const opaque object containing enough information for the underlying 
implementation to determine a unique entity.  This is used to distinquish
which subroutines can change the opaque object and which cannot.
 */
  typedef const void * cEntity_Workset_Iterator;


//====================================================================
// Typedef MeshError
//====================================================================

/**
An opaque object containing information about errors generated in a TSTT
method call.  This will be NULL for successful calls; non-null if an error
is generated.  Methods for checking the error status are provided.
 */
  typedef void * MeshError;

//====================================================================
// Mesh_Create
//====================================================================

/**
   Create a TSTT mesh.  This method must be called before any other
   TSTT methods.

   @param mesh (out) the pointer to the TSTT mesh interface to be created.  No
information will be loaded until Mesh_Load is called.
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Mesh_Create(Mesh_Handle *mesh, MeshError *errorReturn);


//====================================================================
// Mesh_Load
//====================================================================

/**
   Load information into a TSTT mesh.  This method must be called after
   Mesh_Create and before any other methods can be called to access mesh
   or entity information.  Note that multiple Mesh_Load calls can be made
   for a given Mesh_Handle.  It is assumed that the union of information
   from each load call will be available in Mesh_Handle in the other 
   TSTT methods.

   @param mesh (in) the handle to the TSTT mesh for which information is
to be loaded.
   @param name (in) a string identifier for which mesh to load, for example
from a file.
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Mesh_Load(cMesh_Handle mesh, const char *name, MeshError *errorReturn);


//====================================================================
// Mesh_Services_GetInt
//====================================================================

/**
   A services method for getting integer information based on a 
string key value.
     String values that are currently supported include
        "number of 0D entities" or "number of vertices"
        "number of 1D entities" or "number of edges"
        "number of 2D entities" or "number of faces"
        "number of 3D entities" or "number of regions"

   It is anticipated that this function can be easily expanded
   to include other required key value pairs. 
 
   @param mesh (in) the handle to the TSTT mesh for which information is
requested.
   @param key (in) the string identifier for the requested information
   @param value (out) the returned value
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Mesh_Services_GetInt(cMesh_Handle mesh, const char* key, 
                    Int *value, MeshError *errorReturn);


//====================================================================
// Mesh_Services_GetFloat
//====================================================================

/**
   A services method for getting float information based on a 
string key value.  No string key values for this function are supported 
at this time.

   @param mesh (in) the handle to the TSTT mesh for which information is
requested.
   @param key (in) the string identifier for the requested information
   @param value (out) the returned value
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Mesh_Services_GetFloat(cMesh_Handle mesh, const char* key, 
                    Float *value, MeshError *errorReturn);


//====================================================================
// Mesh_Services_GetDouble
//====================================================================

/**
   A services method for getting double information based on a 
string key value.  No string key values for this function are supported 
at this time.

   @param mesh (in) the handle to the TSTT mesh for which information is
requested.
   @param key (in) the string identifier for the requested information
   @param value (out) the returned value
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Mesh_Services_GetDouble(cMesh_Handle mesh, const char* key, 
                    Double *value, MeshError *errorReturn);


//====================================================================
// Mesh_Services_GetString
//====================================================================

/**
   A services method for getting string information based on 
a string key value.  No string key values for this function are supported 
at this time.

   @param mesh (in) the handle to the TSTT mesh for which information is
requested.
   @param key (in) the string identifier for the requested information
   @param value (out) the returned value
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Mesh_Services_GetString(cMesh_Handle mesh, const char* key, 
                    char **value, MeshError *errorReturn);


//====================================================================
// Mesh_GetGeometricDimension
//====================================================================

/**
  Gets the geometric dimension of the mesh.  Note that this does
  not necessarily equal the highest topological entity as surface
  meshes consisting of two-dimensional topological elements can
  live in three space.

   @param mesh (in) the handle to the TSTT mesh for which information is
requested.
   @param dim (out) the integer dimension of the mesh
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/    
   void Mesh_GetGeometricDimension(cMesh_Handle mesh, int *dim, 
                                   MeshError *errorReturn);


//====================================================================
// Mesh_GetEntities
//====================================================================

/**
   Gets an array of entity handles representing all of the entities of type
   entity_type.

   @param mesh (in) the handle to the TSTT mesh for which information is
requested.
   @param entity_type (in) one of the enumerated types given 
             TSTT_VERTEX,  TSTT_EDGE,  TSTT_FACE,  TSTT_REGION
   @param entity_handles (out) the array of opaque entity handles
   @param num_entities (inout) the number of entities in the array
      If num_entities=0 on input, the mesh allocates the entity_handles array.
      If num_entities>=actual on input, it is expected that the application 
        passes in an array of that size, and the mesh fills it.
      If num_entities<actual, an error is returned.
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/    
  void  Mesh_GetEntities(cMesh_Handle mesh, const EntityType entity_type, 
          Entity_Handle **entity_handles, Int *num_entities, 
          MeshError *errorReturn);
    

//====================================================================
// Mesh_GetEntityAdjacency
//====================================================================

/**
   Get the adjacenty list for the entities of the TSTT Mesh.  The entities
   are returned in the form of an array of EntityHandles which uniquely 
   identify each mesh entity of a given type.  We return the data in 
   compressed sparse row (CSR) format.  In this format three arrays are
   returned:
      csr_index: an integer array giving the location in the index array
                 beginning each entity's adjacency list.  For example, 
                 the starting point of the adjacencies for entity k (as
                 given in csr_adjacencies array) is given by csr_index[k], 
                 and the number of adjacencies is given by 
                 csr_index[k+1]-csr_index[k].
      csr_adjacencies: an integer array giving the locations of the adjacent
                 entities in the entity_handles array 
      entity_handles: the array containing unique identifiers for the 
                 requested adjacent entities

   For example to request the vertex adjacency information of two triangles
   and a quadrilateral (the two triangles share vertices 1 and 2, the 
   quadrilateral shares vertices 2 and three whith the second triangle)
      The csr_index array is 0, 3, 6, 10
      The csr_adjacencies array is 0 1 2 1 2 3 2 3 4 5
      The entity_handles array contains 6 elements uniquely identifying the
          vertices of the elements

   @param mesh (in) the handle to the TSTT mesh for which information is
requested.
   @param entity_requestor (in) the type of entity making the adjacency 
            request; one of VERTEX, EDGE, FACE, or REGION
   @param adjacency_entity_requested (in) the entity type requested; one of
            VERTEX, EDGE, FACE, or REGION
   @param entity_handles (out) the array of unique entity handles 
             for the requested adjacencies
   @param csr_pointer (out) an integer array that gives the number of entities
            adjacent to each requestor.  The TSTT mesh is required to allocate 
            this array.
   @param csr_data (out) an integer array giving the locations of 
            the adjacent entities in the entity_handles array.  The TSTT 
            mesh is required to allocate this array.
   @param num_adj_entities (inout) the number of entities in the array
      If num_adj_entities=0 on input, the mesh allocates the entity_handles array.
      If num_adj_entities>=actual on input, it is expected that the application 
        passes in an array of that size, and the mesh fills it.
      If num_adj_entities<actual, an error is returned.
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.

  The following rules apply for what is returned from this function.

  Downward adjacencies are clear.

  Equal adjacencies are not well defined and return an error.  We will
  eventually define neighbor relationships that will fill this gap.

  Upward adjacencies are clear for conformal meshes, non conformal
  meshes are not defined at this time.

  In particular, for the following entity_requestor:adjacency_entity_requested
    pairs we have

    REGION:VERTEX returns the corner vertices of an region

    REGION:EDGE returns the edges defining the region 

    REGION:FACE returns the faces defining the region 

    FACE:VERTEX returns the corner vertices of the face

    FACE:EDGE returns the edges defining the face
 
    EDGE:VERTEX returns the endpoints of the edge

    VERTEX:REGION returns the regions for which it is a corner vertex

    VERTEX:FACE returns the faces for which it is a corner vertex

    VERTEX:EDGE returns the edges for which it is an endpoint

    EDGE:REGION returns the regions which it helps define in that
          its endpoints are corner vertices of the returned regions

    EDGE:FACE returns the faces which it helps define in that its
          endpoints are the corner vertices of the returned faces
   
    FACE:REGION returns the elements which it helps define in that
          its corner vertices are the corner vertices of the returned
          elements 


*/
  void  Mesh_GetEntityAdjacency(cMesh_Handle mesh, 
                const EntityType entity_requestor,
                const EntityType adjacency_entity_requested,
                Entity_Handle **entity_handles, 
                Int **csr_pointer, Int **csr_data,
                Int *num_adj_entities,
                MeshError *errorReturn);


//====================================================================
// Mesh_GetNativeStorageOrder
//====================================================================

/**   
   Gets the preferred storage ordering for the TSTT mesh component. 

   @param mesh (in) the handle to the TSTT mesh for which information is
        requested.
   @param storage_order (out) the preferred storage order for the TSTT
        mesh component.  If no storage order is preferred UNDETERMINED 
        is returned.
   @param errorReturn (out) the error object; will be NULL if successful,
        non-NULL otherwise.
*/
  void Mesh_GetNativeStorageOrder(cMesh_Handle mesh, StorageOrder *storage_order,
                         MeshError *errorReturn);


//====================================================================
// Mesh_GetCoordinates
//====================================================================

/**   
   Get the coordinates of the vertices.

   Return an array of vertex coordinates that will be indexed with a
   global or local ID. The user may choose the storage order to be
   either interleaved, xyzxyzxyz storage or blocked, xxxyyyzzz
   storage.  It is assumed that the coordinates are returned in a
   fixed order that will be reproduced if no mesh modifications have
   occurred between calls.

   @param mesh (in) the handle to the TSTT mesh for which information is
        requested.

   @param num_entities (in/out) the number of entities for which the user is
        requesting coordinate information.  If this argument is zero on input,
        the TSTT mesh component will allocate the coordinate array and fill
        it.  If this argument is not zero, it is assumed that the array
        coordinates is of sufficient length to store the requested coordinate
        information.  If the array is not big enough an error is issued.

   @param coordinates (input/output): If this array is allocated on 
         input then the results will be copied to this array if it is 
         big enough, otherwise an error is issued. If this array 
         has not been allocated space on input then the space will 
         be allocated by this routine.

   @param storage_order (input/output): If this is set to BLOCKED or INTERLEAVED
         on input, the TSTT mesh is required to return the coordinates using
         that storage order.  If this argument is set to UNDETERMINED, the TSTT
         mesh is free to choose one method or the other and set the value of
         storage_order to be the chosen type.
      
   @param options (input) : string containing options separated 
       by semi-colons, options="option1; option2; ..."
         option = "interleaved" : vertices will be returned as xyzxyzxyz...
         option = "blocked"     : vertices will be returned as 
                                  xxxx..xyyyy...yzzzz...z
         option = "range(n1:n2)" : return the vertices in the range 
                   n1,n1+1,...,n2
         option = "local"       : return the vertices local to this processor
         option = "includeHigherOrderNodes" : return nodes used by higher 
                   order elements -- this option needs to be worked out 
                   more carefully to specify exactly where the higher-order 
                   nodes sit etc.
    Example: options="interleaved; range(0:1000)" will return vertices 
             0 to 1000 in the interleaved format.

   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void Mesh_GetCoordinates(cMesh_Handle mesh, Int* num_entities,
                         Double** coordinates, StorageOrder *storage_order,
                         char *options, MeshError *errorReturn);


//====================================================================
// Mesh_GetEntityIndicies
//====================================================================

/**
Return the indicies of the vertices that define all entities of 
a given type.
 
The index values (i.e. the global or local ID of a vertex) for 
       entity i are given by (assuming a compressed-row format) :
            index[e1],...,index[e2-1] 
        where 
            e1=count[i], e2=count[i+1]
   
  @param mesh (in) the handle to the TSTT mesh for which information is
requested.
  @param entityType (in) one of the enumerated types 
          VERTEX, EDGE, FACE, REGION.  This is the entity which is requesting
          the vertex adjacency information
  @param count (input): array to place results in. Usual behaviour 
         depending on whether array has been allocated or not (see above).
  @param index (input) : array to place results in. Usual behaviour 
         depending on whether array has been allocated or not (see above).
  @param options (input) : a list of options, separated by semi-colons 
         options="option1; option2; ...".
         option = "csr" : return in compressed-row format
         option = "expanded" : return in uncompressed format
         option = "HEXAHEDRON" : return only those regions that are 
                   also hexahedra.
         option = "includeHigherOrderNodes" : return nodes used by higher 
                   order elements -- this option needs to be worked out 
                   more carefully to specify exactly where the 
                   higher-order nodes sit etc.
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void Mesh_GetEntityIndicies(cMesh_Handle mesh, 
                   const EntityType entityType, 
                   Int** count, Int* num_count, 
                   Int** index,  Int* num_indicies,
                   char *options, MeshError *errorReturn );


//====================================================================
// Mesh_GetEntityAdjacencyIndicies
//====================================================================

/**
Return adjacency information.

  Example: To return the list of faces that are adjacent to all regions 
     one would set entityType=REGION, entityAdjacencyType=FACE.

     The index values for the adjacency i (numbered according to the 
     order they are returned in the call to getEntitiesIndicies) are 
     given by (assumming a compressed-row format)
            index[e1],...,index[e2-1] 
        where 
            e1=count[i], e2=count[i+1]

  @param mesh (in) the handle to the TSTT mesh for which information is
requested.
  @param entityType (in): the requesting entity type, one of VERTEX,
        EDGE, FACE, or REGION
  @param entityAdjacencyType (in): the requested entity type adjacency, 
        one of VERTEX, EDGE, FACE, or REGION.  See the documentation for
        Mesh_GetEntityAdjacency for the rules describing the adjacency
        behavior.
  @param count (inout):  array to place results in. Usual 
        behaviour depending on whether array has been allocated or not 
        (see above).
  @param index (out): array to place results in. Usual behaviour 
        depending on whether array has been allocated or not (see above).
  @param options (in) : a list of options, separated by semi-colons, 
        options="option1; option2; ...".
         option = "csr" : return in compressed-row format
         option = "expanded" : return in uncompressed format
         option = "includeHigherOrderNodes" : return nodes used by 
                   higher order elements -- this option needs to be 
                   worked out more carefully to specify exactly where 
                   the higher-order nodes sit etc.
 @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void Mesh_GetEntityAdjacencyIndicies(cMesh_Handle mesh, 
               const EntityType entityType, 
               const EntityType entityAdjacencyType,
               Int **count, Int* num_count, 
               Int **index, Int* num_indices,
               char *options, MeshError *errorReturn);



//====================================================================
// Mesh_InitializeWorksetIterator
//====================================================================

/** 
   Initialize the entity workset iterator.  The entities are returned
   in the Mesh_GetNextWorkset function defined below.  Each iterator
   has an entity type and preferred workset size associated with it.
   Note that the Mesh_GetEntities function can be replaced with these
   two calls if desired.

   @param mesh (in) the handle to the TSTT mesh for which the iterator
     will be defined
   @param entity_type (in) one of the enumerated types given 
         VERTEX, EDGE, FACE, REGION
   @param workset_size (in) the number of entities in each workset.  
         The entities are returned in chunks of this size until the
         last workset for which the actual size will be less than
         requested size.  Note that this can be set to a default of
         -1 to get the entire list of entities of a particular type.
   @param workset_iterator (out) the workset iterator.  This structure
         should keep enough state information that the function 
         Mesh_GetNextWorkset will be efficient.  Note, it is assumed
         that workset_iterator keeps a pointer to its parent mesh
         object so that once it's created it cannot be used with 
         multiple meshes.
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
   void Mesh_InitializeWorksetIterator(cMesh_Handle mesh, 
	      const EntityType entity_type, 
              const Int requested_workset_size,
              Entity_Workset_Iterator *workset_iterator,
              MeshError *errorReturn);

//====================================================================
// Mesh_GetNextWorkset
//====================================================================

/** Iterator to get the next workset.  Workset_size entities are 
    returned in the entity handle array until the final workset. 
    This allows us to traverse the entities in manageable chunks
    and can be tuned to provide a balance between memory overhead
    and performance of accessing information through interfaces.

   @param workset_iterator (in) the workset iterator.  This structure
         should keep enough state information that entity access
         provided by this function is efficient. Note, it is assumed
         that workset_iterator keeps a pointer to its parent mesh
         object so that a mesh object does not need to be passed in
         to this function.
   @param entity_handles (out) the array of opaque entity handles
         in the workset
   @param num_entities (out) the number of entities in the workset
         array. num_entities will usually equal workset_size until 
         the last workset is reached in which case it will be less 
         than or equal to worksetsize.  If the last workset contains
         a full worksetsize number of entities, an additional call to
         this function will return 0 entities.
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
    void Mesh_GetNextWorkset(cMesh_Handle mesh, 
           Entity_Workset_Iterator workset_iterator,
           Entity_Handle **entity_handles, Int *num_entities, 
           MeshError *errorReturn);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  "WorkSet" is not defined in AOMD, as a temporary solution, we provide
//   the following two functions, which are equivalent to the above   
//   functions but do the job at entity level.
//====================================================================
// Mesh_InitializeEntityIterator (proposed and iimplemented by yluo)
//====================================================================
/** This function initializes the iterator for "entityType". 
    @param mesh (in): the mesh on which the iterator to be initialized
    @param entityType (in): the entity type the iterator to iterate
    @param firstEntity (out): the first entity of type "entityType"
    @param lastEntity  (out): the last entity of type "entityType"
    @param errorReturn (out): error message 
 */
   void Mesh_InitializeEntityIterator(cMesh_Handle mesh, 
               EntityType entityType, 
	       Entity_Handle *firstEntity,
	       Entity_Handle *lastEntity, 
	       MeshError *errorReturn);
	   
//====================================================================
// Mesh_EntityGetNextIterator (proposed and implemented by yluo)
//====================================================================
/** move the iterator, whcih points to one entity of type "entityType", 
    to point to the next entity of type "entityType".
    @param mesh (in): the mesh the iterator is iterating.
    @param entityType (in): the entity type the iterator is iterating.
    @param entity (in/out): an iterator which is pointing to an 
           entity of type "entityType".
    @param errorReturn (out): error message.
 */
   void Mesh_GetNextEntity(cMesh_Handle mesh, EntityType entityType, 
           Entity_Handle *entity, MeshError *errorReturn);

//====================================================================
// Mesh_FreeEntityHandles
//====================================================================

/**
  Frees the information in the entity_handles array.  Note that this
does not delete the entities themselves, merely the opaque pointers
to them.
  @param entity_handles (in) the array to free. 
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
 */    
  void  Mesh_FreeEntityHandles(cMesh_Handle mesh, 
                Entity_Handle *entity_handles, 
                MeshError *errorReturn);

//====================================================================
// Mesh_FreeEntityAdjacency
//====================================================================

/**
  Frees the information in the entity_adjacencies arrays.  Note that this
does not delete the adjacent entities, merely the index arrays and
opaque pointers associated with them.

  @param entity_handles (in) the array of opaque entity handles to free.  
  @param csr_pointer (in) an integer array giving the number of entities
            adjacent to each requestor to free. 
  @param csr_data (in) an integer array giving the locations of 
            the adjacent entities in the entity_handles array to free. 
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
 */
  void  Mesh_FreeEntityAdjacency(cMesh_Handle mesh,
                Entity_Handle *entity_handles, 
                Int *csr_pointer, Int *csr_data, 
                MeshError *errorReturn);


//====================================================================
// Mesh_Destroy
//====================================================================

/** Destroys the mesh information that is referenced by the opaque mesh
handle mesh.

  @param mesh (in) the handle to the TSTT mesh for which information is
requested.
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
 */
  void  Mesh_Destroy(Mesh_Handle mesh, MeshError *errorReturn);

//====================================================================
// Mesh_GetEntitiesFromID
//====================================================================

/** OPTIONAL FUNCTION that returns an array of opaque entity handles from an
array of unique global IDs.
  @parem mesh (in) the mesh from which the entity will come
  @param ID (in) an array of unique entity identifiers
  @param entity_handles (out) an array of the corresponding opaque entity handles,
       it is assumed that the TSTT mesh allocates this array
  @param errorReturn (out) the error object; will be NULL if successful,
       non-NULL otherwise.
*/
  void  Mesh_GetEntitiesFromID(cMesh_Handle mesh, Int *ID, EntityType entity_type, 
                const Int num_identifiers, Entity_Handle **entity_handles, 
                MeshError *errorReturn);

//====================================================================
// Entity_GetGlobalID
//====================================================================

/** OPTIONAL FUNCTION that returns an array of Identifiers for the
input opaque entity_handles.  The Identifier is unique among entites
of the same dimension (should this be among all entities?).
  @param entity_handles (in) an array of opaque entity handles
  @param num_entities (in) the number of entities in the array
  @param ID (out) an array containing the integer IDs; it is assumed
         that the underlying component will allocate the space for the
         Identifier array (does it make sense provide the size of
         Identifier to the application so that the app can manage the
         memory?)
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Entity_GetGlobalID(cMesh_Handle mesh, 
                cEntity_Handle *entity_handles, 
                Int num_entities,
                Int **ID, MeshError *errorReturn);


//====================================================================
// Entity_GetGeometricDimension
//====================================================================

/** Returns an integer array of spatial dimensions for the input array
    of Entity_Handles.  Note that the geometric dimension can be greater 
    than or equal to the topological dimension (for example, vertices or 
    faces embedded in three-dimensional space).
  @param entity_handle (in) an array of opaque entity handles
  @param num_entities (in) the number of entities in the array
  @param dim (out) an array containing the integer spatial 
       dimensions of the entities; it is assumed that the array
       has been allocated by the application to be of size num_entities.
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  /* void  Entity_GetGeometricDimension(cEntity_Handle *entity_handles, 
                const Int num_entities,
                Int *dim, MeshError *errorReturn);*/
    

//====================================================================
// Entity_GetTopologicalDimension
//====================================================================

/** Returns an integer array of topological dimensions for an input
  array of opaque Entity_Handles.
  @param entity_handle (in) an array of opaque entity handles
  @param num_entities (in) the number of entities in the array
  @param dim (out) an array containing the integer topological
       dimensions of the entities; it is assumed that the array
       has been allocated by the application to be of size num_entities.
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Entity_GetTopologicalDimension(cMesh_Handle mesh,
                cEntity_Handle *entity_handles, 
                const Int num_entities,
                Int *dim, MeshError *errorReturn);
    
//====================================================================
// Entity_GetTopology
//====================================================================

/**
  Get the topologies of an array of entities.  The value is 
   returned as a value from the enumerated type EntityTopolgy.  

  @param entity_handles (in) an array of opaque entity handles
  @param num_entities (in) the number of entities in the array
  @param topology (out) an array containing the topology of the entity; it 
       is assumed that the array has been allocated by the application to 
       be of size num_entities.
  @param errorReturn (out) the error object; will be NULL if successful,
       non-NULL otherwise.
*/
  void  Entity_GetTopology(cMesh_Handle mesh,
                cEntity_Handle *entity_handles, 
                const Int num_entities,
                EntityTopology *topology,
                MeshError *errorReturn);

//====================================================================
// Entity_GetType
//====================================================================

/**
  Get the type of an array of entities.  The value is 
   returned as a value from the enumerated type EntityType.  Note that
   duplicates the functionality in Entity_GetTopologicalDimension, but
   allows the user to work directly with the enumerated type rather than
   an integer.

  @param entity_handles (in) an array of opaque entity handles
  @param num_entities (in) the number of entities in the array
  @param type (out) an array containing the entity types; one of 
       VERTEX, EDGE, FACE or REGION.  It is assumed that the array is
       allocated by the application
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Entity_GetType(cMesh_Handle mesh,
                cEntity_Handle *entity_handles, 
                const Int num_entities,
                EntityType *type,
                MeshError *errorReturn);

//====================================================================
// Entity_GetVertexCoords
//====================================================================

/** Get the coordinates of an array of vertex mesh entities.
    The storage order options are handled as an enumerated type
    because it will be faster and more typesafe to check than a 
    string option

  @param entity_handles (in) an array of opaque entity handles
  @param num_entities (in) the number of entities in the array
  @param storage_order (in) the desired order of storage, blocked or
        interleaved  
  @param coords (inout) a double array of size num_entities*num_dimension 
        containing the vertex coordinates
  @param num_coords (inout) an integer giving the number of doubles
     in the array
     If num_coords=0, the mesh allocates the array 
     If num_coords>=dimension*num_entites, the application passes in an array 
        of that size, and the mesh fills it 
     If num_coords<dimension*num_entites requested, an error is returned 
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Entity_GetVertexCoords(cMesh_Handle mesh, 
               cEntity_Handle *entity_handles, Int num_entities, 
               StorageOrder storage_order, Double **coords, 
               Int *num_coords, MeshError *errorReturn);
   
//====================================================================
// Entity_SetVertexCoords  // NON STANDARD // added by tleurent@mcs.anl.gov
//====================================================================

/** Set the coordinates of an array of vertex mesh entities.
    The storage order options are handled as an enumerated type.

  @param entity_handles (in) an array of opaque entity handles
  @param num_entities (in) the number of entities in the array
  @param storage_order (in) the order of storage, blocked or
        interleaved
  @param num_dimension (in) 2 or 3 ... number of coordinates per vertex
  @param coords (in) a double array of size num_entities*num_dimension 
        containing the vertices coordinates
  @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Entity_SetVertexCoords(cMesh_Handle mesh, 
               cEntity_Handle *entity_handles, Int num_entities, 
               StorageOrder storage_order, Int num_dimension,
               Double *coords, MeshError *errorReturn);

//====================================================================
// Entity_GetAdjacencies
//====================================================================

/** 
Returns the number of adjacent entities of a given dimension and 
an array of unique identifiers for the adjacent entities.

   @param entity_handles (in) an array of opaque entity handles
   @param num_entities (in) the number of entities in the array
   @param entity_type (in) the type of desired adjacencies
   @param adj_entities_handle (out) the array of opaque pointers containing
the adjacent entity information
   @param csr_pointer (out) an integer array that gives the number of entities
            adjacent to each requestor; it is assumed that the array is allocated
            by the TSTT mesh
   @param csr_data (out) an integer array giving the locations of 
            the adjacent entities in the entity_handles array; it is assumed
            that the array is allocated by the TSTT mesh
   @param num_adjacent (out) the number of adjacent entities
     If num_adjacent=0, the mesh allocates the adj_entities_handle array 
     If num_adjacent>=actual, the application is assumed to have passed 
in an array  for adj_entities_handle
        of that size, and the mesh fills it 
     If num_adjacent<actual, an error is returned
   @param errorReturn (out) the error object; will be NULL if successful,
non-NULL otherwise.
*/
  void  Entity_GetAdjacencies(cMesh_Handle mesh, 
               cEntity_Handle *entity_handles,
               Int num_entities,
               const EntityType entity_type, 
               Entity_Handle **adj_entities_handle, 
 	       Int **csr_pointer, Int **csr_data,
               Int *num_adjacent,
               MeshError *errorReturn);


/* ============================================================
  PROPOSALS
  ============================================================ */
    
  /* TAG FUNCTIONALITY TO BE REVISED BY TIM */
  void Mesh_tagGetHandle (cMesh_Handle mesh,
                          const char *tag_name,
                     void** tag_handle,
                     MeshError *errorReturn);
  
  void Mesh_GetTag_Entity(cMesh_Handle mesh,
                          cEntity_Handle this_entity,
                          void* tag_handle,
                          void **tag_value, Int *tag_size,
                          MeshError *errorReturn);
			
// The following set of functions are from SIDL version with minor
// modifications. They allow to add a new tag to all elements in a
// mesh, delete a tag from all elements, change the data for a specific
// element and get tag data from a specific element.

void Mesh_AddTag (cMesh_Handle mesh, const char* tag_name,
		  void* default_value,
		  void* tag_handle,
		  MeshError *errorReturn);
		  
void Mesh_DeleteTag (cMesh_Handle mesh, void* tag_handle,  
		      MeshError *errorReturn);
		      
void Entity_SetTag (cMesh_Handle mesh, void* tag_handle,
                    cEntity_Handle entity_handle,
		    void* tag_value,
		    MeshError *errorReturn);
		      
void Entity_GetTag (cMesh_Handle mesh, void* tag_handle,
                    cEntity_Handle entity_handle,
		    void*& tag_value,
		    MeshError *errorReturn);
  /* ITERATORS TO BE REVISED BY OTTMAR */

  /* HIGHER ORDER NODE ACCESS TO BE REVISED BY TIM/PROPOSED 
     BY OTTMAR */  

  /* SUBSET ACCESS TO VERTEX COORDS/ADJ INFO TO BE PROPOSED BY
     LORI/RAY */

  /* SIDL SPEC TO BE PROPOSED BY KYLE */    

/*==============================================================
  NEEDED
  ============================================================ */

  /* ERROR CODES THAT ARE RETURNED FROM EACH FUNCTION */

  /* A SERVICES OBJECT THAT LETS THE USER KNOW WHAT ADJACENCY INFO IS
     SUPPORTED */

/*========================================================================
  ERROR CODE
  ========================================================================*/

typedef void ErrorCallbackFunction(MeshError *error);

enum DefaultActions
{
    silent=0,
    printErrors,
    abortOnErrors
};
  
/* set an error integer code and a description */
  void MeshError_Set(const int err, const char *description);

  void MeshError_SetDefaultAction(MeshError *error, const DefaultActions action );

  void MeshError_GetDefaultAction(MeshError *error, DefaultActions *action);
 
  /* return 1 if there is no error */
  int MeshError_OK(MeshError *error);

  void MeshError_Get(MeshError *error, int *err, char *description);
  int  MeshError_GetErrorNumber(MeshError *error);
  void MeshError_GetErrorDescription(MeshError *error, char *description);
  void MeshError_PrintError(MeshError *error, const char* label);

  /* set the call back function */
  void Mesh_SetErrorCallbackFunction(ErrorCallbackFunction function);

  /* Set the default action to take on errors and warnings */
  void Mesh_SetDefaultAction(DefaultActions action);
}

#endif
