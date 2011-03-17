#ifndef _ITAPS_iMesh
#define _ITAPS_iMesh

  /**\brief iMesh Interface Specification Major Release Number
    */
#define IMESH_MAJOR_VERSION 1

  /**\brief iMesh Interface Specification Minor Release Number
    */
#define IMESH_MINOR_VERSION 1

  /**\brief iMesh Interface Specification Patch Release Number
    *
    * Technically speaking, there is not much practical value in
    * a 'patch' digit for an interface specification. A 'patch'
    * release is typically only used for bug fix releases. That
    * is unlikely to apply to an interface specification such as
    * iMesh. Nonetheless, there is non-technical value in having
    * a bit finer grained control over version numbers that a patch
    * digit affords. So, we define one for iMesh.
    */
#define IMESH_PATCH_VERSION 0

  /** \mainpage The ITAPS Mesh Interface iMesh
   *
   * The ITAPS Mesh Interface iMesh provides a common interface for
   * accessing mesh and data associated with a mesh.  Applications written
   * to use this interface can use a variety of implementations, choosing
   * the one that best meets its needs.  They can also use tools written
   * to this interface, for example mesh smoothing, adaptive mesh refinement,
   * and parallel mesh support.
   *
   * \section ITAPS Data Model
   *
   * The ITAPS interfaces use a data model composed of four basic data types: \n
   * \em Entity: basic topological entities in a mesh, e.g. vertices, 
   * triangles, hexahedra. \n
   * \em Entity \em Set: arbitrary grouping of other entities and sets. 
   * Entity sets also support parent/child relations with other sets which
   * are distinct from entities contained in those sets.  Parent/child links
   * can be used to embed graph relationships between sets, e.g. to 
   * represent topological relationships between the sets. \n
   * \em Interface: the object with which mesh is associated and on which
   * functions in iMesh are called. \n
   * \em Tag: application data associated with objects of any of the other 
   * data types.  Each tag has a designated name, size, and data type.
   *
   * \section ITAPS Entity Type, Topology
   * Each entity has a specific Entity Type and Entity Topology.  The Entity 
   * Type is one of VERTEX, EDGE, FACE, and REGION, and is synonymous with
   * the topological dimension of the entity.  The Entity Topology denotes
   * the specific shape, for example TRIANGLE, QUADRILATERAL, TETRAHEDRON,
   * and HEXAHEDRON.  Entity Type and Entity Topology exist as enumerated
   * types, Entity Type in the iBase_EntityType enumeration, and
   * Entity Topology in the iMesh_EntityTopology enumeration.
   *
   * \section ITAPS Entity-, Array-, and Iterator-Based Access
   *
   * The iMesh interface provides functions for accessing entities
   * individually, as arrays of entities, or using iterators.  These access
   * methods have different memory versus execution time tradeoffs, 
   * depending on the implementation.
   *
   * \section ITAPS Lists Passed Through Interface
   *
   * Many of the functions in iMesh have arguments corresponding to lists of 
   * objects.  In-type arguments for lists consist of a pointer to an array and
   * a list size.  Lists returned from functions are passed in three arguments,
   * a pointer to the array representing the list, and pointers to the
   * allocated and occupied lengths of the array.  These three arguments are 
   * inout-type arguments, because they can be allocated by the application and
   * passed into the interface to hold the results of the function.  Lists
   * which are pre-allocated must be large enough to hold the results of the
   * function; if this is not the case, an error is generated.  Otherwise, the
   * occupied size is changed to the size output from the function.  If a list
   * argument is unallocated (the list pointer points to a NULL value) or if
   * the incoming value of the allocated size is zero, the list storage will be
   * allocated by the implementation.  IN ALL CASES, MEMORY ALLOCATED BY ITAPS
   * INTERFACE IMPLEMENTATIONS IS DONE USING THE C MALLOC FUNCTION, AND CAN BE
   * DE-ALLOCATED USING THE C FREE FUNCTION.
   *
   * \section ITAPS Storage Orders
   *
   * Many of the functions in iMesh can return arrays of tuples; that is, arrays
   * of multi-valued type. For example, the function iMesh_getVtxArrCoords,
   * returns an array of xyz coordinate 3-tuples (or, perhaps for geometrically
   * 2D meshes, xy 2-tuples). In these situations, there are multiple ways the
   * data can be organized in memory. For example, it could be stored xyz,xyz,xyz
   * or xxx...,yyy...,zzz.... These two different storage orders are referred
   * to as INTERLEAVED and BLOCKED, respectively. For some functions in iMesh,
   * the storage order is explicitly specified as an argument to the function.
   * For other functions, the storage order is not explicitly specified. And,
   * in these cases, it shall always be implicitly assumed to be INTERLEAVED.
   * This fact will be mentioned in the documentation for each specific function
   * where it applies. For example, in the case of iMesh_getEntArrAdj, the returned
   * array of adjacent entities is multi-valued in that it stores for each
   * entity queried, all its adjacent entities. Such an array will be stored
   * INTERLEAVED with all adjacent entities for the first entity in the query
   * followed by all adjacent entities for the second entity in the query and
   * so forth.
   *
   * \section String length arguments
   *
   * Many of the functions in iMesh involve passing a string and also the length
   * of that string. The question arises, how is the null character ('\0') handled?
   * For users of the iMesh interface calling iMesh functions, it is optional
   * as to whether or not to include the nil character in computing the length of
   * the string. So, for example, calling iMesh from a C program, users could
   * pass strlen(my_string) or strlen(my_string)+1 as the length of the string.
   *
   * However, it should be noted that the situation is different for implementers
   * of the iMesh interface. In implementing an iMesh interface function, there
   * can be no assumption that the string is indeed null terminated. The length
   * argument the caller passes in may or may NOT include the null character and
   * implementations must be coded to accommodate this. This requirement is 
   * primarily due to differences in how Fortran and C/C++ handle passing of
   * strings as function arguments.
   *
   * \section Cycles in Set-Inclusion and Parent-Child structures.
   *
   * There are two graph-like structures in the iMesh interface and data
   * model; the set-inclusion structure and the parent-child link structure.
   * Whether these structures support cycles is relevant to implementors.
   *
   * Over the evolution of the iMesh data model and API, both of these
   * structures have been viewed more or less like a tree and so cycles seem
   * incompatible with that notion.
   *
   * Allowing a cycle in the set inclusion structure implies all entity sets
   * in the cycle are all equal to each other. That is the only rational,
   * model-level view that would allow them all to be (improper) subsets of
   * each other. On the other hand if the iMesh specification excludes cycles
   * from the set inclusion structure, the time complexity (performance) as a
   * function of the number of entity sets may be prohibitive for
   * implementations to detect and prevent them.
   *
   * Allowing a cycle in the parent-child link structure seems likewise hard
   * to justify. However, when the parent-child structure is viewed more like
   * a general graph (a view that the current API itself supports even if the
   * function names themselves do not suggest that) than specifically a tree,
   * the admission of cycles there is potentially more natural and useful.
   *
   * Implementations are required to support cycles in the Parent-Child
   * structure. Implementations are neither required to support nor required
   * to explicitly prevent cycles in the Set-Inclusion structure. Portable
   * applications should NOT rely on implementations support for cycles
   * in the set-inclusion structure.
   *
   * \section Order-preserving versus duplicate-preventing entity sets
   *
   * An entity set in iMesh supports managing its contained members in one of
   * two modes. When an entitiy set is first created, the caller is required
   * to indicate which mode of membership the set will support. The two modes
   * that are supported are a) order-preserving or b) duplicate-preventing.
   *
   * For order-preserving membership, the implementation will permit duplicate
   * members. However, the implementation will guarantee that the order in which
   * members are added to the set will be the same as the order in which they are
   * queried by the various methods that return the members of a set. This order
   * preserving guarantee holds across removals. However, it does not hold across
   * removals followed by re-additions of the previously removed members. This
   * kind of an entity set behaves like an STL vector or STL list.
   *
   * For duplicate-preventing membership, the implementation will guarantee that
   * duplicate members are prevented. Any attempts to add duplicate members to
   * such a set will be detected, prevented and silently ignored by the
   * implementation.  This kind of entity set behaves like an STL set.
   *
   * Finally, although entity sets may contain two kinds of members, entities
   * or other entity sets, the above statements hold only for the entity members.
   * Behavior for entity set members is presently implementation defined.
   *
   * \section Indirection in queries on Set-Inclusion and Parent-Child structures.
   *
   * Various functions to query entities, entity sets and parent or child sets 
   * as well as the numbers of these involve a num_hops argument. If the set
   * upon which the query is originated is the root set, the num_hops argument
   * is irrelevant and is ignored by the implementation. Otherwise, the num_hops
   * argument represents the maximum number of levels of indirection employed in
   * satisfying the query not including the originating set. For example, using
   * value for num_hops of 0 (zero) in iMesh_getEntSets will return all the 
   * entity sets that are immediately contained in a given set. Likewise, a
   * value for num_hops of 1 (one) will return all entity sets that are 
   * immediately contained in the given set plus all entity sets that
   * are contained in those immediately contained sets (e.g. one level of
   * indirection). Using a value of -1 for num_hops will return results for
   * all possible levels of indirection. In other words, using a value of
   * -1 for num_hops is equivalent to setting the maximum number of levels
   * of indirection to infinity.
   *
   * \section Options in newMesh, save and load methods.
   *
   * A few of the functions in iMesh support arbitrary options passed as a
   * character string, called an 'Option String'. The format of and handling
   * of an Option String is as follows...
   *    - Option Strings are INsensitive to case.
   *    - Each option in an Option String is pre-pended with the implementation
   *      name followed by a special character called the separator character.
   *    - The separator is a colon, ':'.
   *    - Multiple options existing in a single Option String are separated by a
   *      special character called the delimiter character.
   *    - The delimiter character is a space, ' '.
   *    - The effect of multiple options in a single Option String is 
   *      INsensitive to order of occurrence in the string.
   *    - By default, implementations silently ignore any options that
   *      do not match on the implementation name part (everything before
   *      the separator character). This way, a caller may included options
   *      in a single string intended for multiple different implementations.
   *    - Implementations may (or may not) warn or error for option strings
   *      that match on implementation name part but are found to be in error
   *      for other reasons the implementation decides.
   *    - Whenever either the separator character, ':', or delimiter character,
   *      ' ', need to appear in an option, they must be escaped with the
   *      backslash character, '\'.
   *
   *  For example, consider the Options String
   *      "grummp:silant FMDB:TwoPhaseIO moab:mpiio_hints\ foo\:bar"
   *
   *  In the above example, the space serves as the delimiter character
   *  between multiple options in the string. The colon serves as the
   *  implementation-name/option separator character. Because options are
   *  required to be insensitive to case, the caller is free to use case as a
   *  word separator as in 'TwoPhaseIO' and even in the implementation name,
   *  as in 'FMDB:', although 'fmdb:twophaseio' and 'fmdb:TWOPHASEIO' would
   *  all have the same effect. In the moab option, both the separator
   *  character and delimiter character appear in the option and so are
   *  pre-pended (e.g. escaped) with the backslash character.

   *  GRUMMP will silently ignore the FMDB: and moab: options because they do
   *  NOT match on the implementation name part. However, GRUMMP may
   *  optionally error out, or warn or silently ignore 'grummp:silant' (it was
   *  supposed to be spelled 'silent') as an invalid option.
   *
   *  Note that iMesh itself currently does not define any options. In order
   *  to discover options a given implementation defines, users are directed
   *  to the developers of the respective implementations.
   */

#include "iBase.h"
#include "iMesh_protos.h"

/**\brief Useful macro for greater-than or equal-to comparisons of version number.
 *
 * Returns true (non-zero) if the current version is greater-than or equal-to that
 * specified by Maj,Min,Pat triple.
 */
#define IMESH_VERSION_GE(Maj,Min,Pat) \
    (((IMESH_MAJOR_VERSION==(Maj)) && (IMESH_MINOR_VERSION==(Min)) && (IMESH_PATCH_VERSION>=(Pat))) || \
     ((IMESH_MAJOR_VERSION==(Maj)) && (IMESH_MINOR_VERSION>(Min))) || \
      (IMESH_MAJOR_VERSION>(Maj)))

#define IMESH_VERSION_TAG__(X,Y,Z) iMesh_Version_##X##_##Y##_##Z
#define IMESH_VERSION_TAG_(X,Y,Z) IMESH_VERSION_TAG__(X,Y,Z)
#define IMESH_VERSION_TAG IMESH_VERSION_TAG_(IMESH_MAJOR_VERSION,IMESH_MINOR_VERSION,IMESH_PATCH_VERSION)

#define IMESH_VERSION_STRING__(X,Y,Z) "iMesh_Version_" #X "." #Y "." #Z
#define IMESH_VERSION_STRING_(X,Y,Z) IMESH_VERSION_STRING__(X,Y,Z)
#define IMESH_VERSION_STRING IMESH_VERSION_STRING_(IMESH_MAJOR_VERSION,IMESH_MINOR_VERSION,IMESH_PATCH_VERSION)

#define IMESH_NEW_MESH_NAME__(A,B,C) A##_##B##_##C
#define IMESH_NEW_MESH_NAME_(A,B,C) IMESH_NEW_MESH_NAME__(A,B,C)
#define IMESH_NEW_MESH_NAME(A) IMESH_NEW_MESH_NAME_(A,IMESH_MAJOR_VERSION,IMESH_MINOR_VERSION)

/*
#undef  iMesh_newMesh
#define iMesh_newMesh IMESH_NEW_MESH_NAME(iMesh_newMesh)
*/

#ifdef __cplusplus
extern "C" {
#endif

    /**\brief  Type used to store iMesh interface handle
     *
     * Type used to store iMesh interface handle
     */
  typedef struct iMesh_Instance_Private* iMesh_Instance;

    /**\brief  Enumerator specifying entity topology
     *
     * Enumerator specifying entity topology.
     */
  enum iMesh_EntityTopology {
    iMesh_EntityTopology_MIN = 0,
    iMesh_POINT = iMesh_EntityTopology_MIN, /**< a general zero-dimensional entity  */
    iMesh_LINE_SEGMENT,       /**< a general one-dimensional entity  */
    iMesh_POLYGON,            /**< a general two-dimensional element  */
    iMesh_TRIANGLE,           /**< a three-sided, two-dimensional element  */
    iMesh_QUADRILATERAL,      /**< a four-sided, two-dimensional element  */
    iMesh_POLYHEDRON,         /**< a general three-dimensional element */
    iMesh_TETRAHEDRON,        /**< a four-sided, three-dimensional element whose
			       *     faces are triangles */
    iMesh_HEXAHEDRON,         /**< a six-sided, three-dimensional element whose
			       *     faces are quadrilaterals */
    iMesh_PRISM,              /**< a five-sided, three-dimensional element which
			       *     has three quadrilateral faces and two
			       *     triangular faces  */
    iMesh_PYRAMID,            /**< a five-sided, three-dimensional element
			       *     which has one quadrilateral face and four
			       *     triangular faces */
    iMesh_SEPTAHEDRON,        /**< a hexahedral entity with one collapsed edge */
    iMesh_ALL_TOPOLOGIES,     /**< allows the user to request information
			       *     about all the topology types */
    iMesh_EntityTopology_MAX = iMesh_ALL_TOPOLOGIES
  };

    /**\brief  Get the error type returned from the last iMesh function
     *
     * Get the error type returned from the last iMesh function.  Value
     * returned is a member of the iBase_ErrorType enumeration.
     * \param instance iMesh instance handle
     * \param *error_type Error type returned from last iMesh function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getErrorType(iMesh_Instance instance,
                          /*out*/ int *error_type, 
                          int *err);

    /**\brief  Get a description of the error returned from the last iMesh function
     *
     * Get a description of the error returned from the last iMesh function
     * \param instance iMesh instance handle
     * \param descr Pointer to a character string to be filled with a
     *        description of the error from the last iMesh function
     * \param *err Pointer to error type returned from function
     * \param descr_len Length of the character string pointed to by descr
     */
  void iMesh_getDescription(iMesh_Instance instance,
                            /*inout*/ char *descr, 
                            int *err, 
                            /*in*/ int descr_len);

    /**\brief  Construct a new iMesh instance
     *
     * Construct a new iMesh instance, using implementation-specific
     * options
     * \param options Pointer to implementation-specific options string
     * \param instance Pointer to iMesh instance handle returned from function
     * \param *err Pointer to error type returned from function
     * \param options_len Length of the character string pointed to by options
     */
  void iMesh_newMesh(const char *options,
                     /*out*/ iMesh_Instance *instance, 
                     /*out*/ int *err, 
                     /*in*/ int options_len);

    /**\brief  Destroy an iMesh instance
     *
     * Destroy an iMesh instance
     * \param instance iMesh instance to be destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_dtor(iMesh_Instance instance, 
                  /*out*/ int *err);

    /**\brief  Load a mesh from a file
     *
     * Load a mesh from a file.  If entity set is specified, loaded mesh
     * is added to that set; specify root set if that is not desired.
     * \param instance iMesh instance handle
     * \param entity_set_handle Set to which loaded mesh will be added, 
     *        root set if not desired
     * \param name File name from which mesh is to be loaded
     * \param options Pointer to implementation-specific options string
     * \param *err Pointer to error type returned from function
     * \param name_len Length of the file name character string
     * \param options_len Length of the options character string
     */
  void iMesh_load(iMesh_Instance instance,
                  /*in*/ const iBase_EntitySetHandle entity_set_handle,
                  /*in*/ const char *name, 
                  /*in*/ const char *options,
                  /*out*/ int *err, 
                  /*in*/ int name_len, 
                  /*in*/ int options_len);

    /**\brief  Save a mesh to a file
     *
     * Save a mesh to a file.  If entity set is specified, save only the
     * mesh contained in that set.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being saved
     * \param name File name to which mesh is to be saved
     * \param options Pointer to implementation-specific options string
     * \param *err Pointer to error type returned from function
     * \param name_len Length of the file name character string
     * \param options_len Length of the options character string
     */
  void iMesh_save(iMesh_Instance instance,
                  /*in*/ const iBase_EntitySetHandle entity_set_handle,
                  /*in*/ const char *name, 
                  /*in*/ const char *options,
                  /*out*/ int *err, 
                  /*in*/ const int name_len, 
                  /*in*/ int options_len);

    /**\brief  Get handle of the root set for this instance
     *
     * Get handle of the root set for this instance.  All mesh in
     * this instance can be accessed from this set.
     * \param instance iMesh instance handle
     * \param root_set Pointer to set handle returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getRootSet(iMesh_Instance instance,
                        /*out*/ iBase_EntitySetHandle *root_set, 
                        /*out*/ int *err);

    /**\brief  Get the geometric dimension of mesh represented in this instance
     *
     * Get the geometric dimension of mesh represented in this instance
     * \param instance iMesh instance handle
     * \param geom_dim Pointer to dimension returned from this function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getGeometricDimension(iMesh_Instance instance,
                                   /*out*/ int *geom_dim, 
                                   /*out*/ int *err);

   /**\brief Set geometric dimension of vertex coordinates
    *
    * Set the geometric dimension of the mesh.  Notes:  An application
    * should not expect this function to succeed unless the mesh database
    * is empty (no vertices created, no files read, etc.)
    *\param instance Mesh database from which to request change.
    *\param geom_dim Requested geometric dimension.
    */
  void iMesh_setGeometricDimension( iMesh_Instance instance,
                                    /*in*/ int geom_dim,
                                    /*out*/ int* err );

    /**\brief  Get the default storage order used by this implementation
     *
     * Get the default storage order used by this implementation.  Value
     * returned is a member of the iBase_StorageOrder enumeration.
     * \param instance iMesh instance handle
     * \param order Pointer to storage order returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getDfltStorage(iMesh_Instance instance,
                            /*out*/ int *order, 
                            /*out*/ int *err);

    /**\brief  Get the adjacency table for this implementation
     *
     * Get the adjacency table for this implementation.  This table 
     * is a 4x4 array whose entries characterize how the implementation
     * behaves when adjacent and intermediate entities are queried.
     * Entry [i,j] (i=row, j=col, 0-based indexing) will have one of
     * the values in the iBase_AdjacencyCost enumeration. Off-diagonal
     * entres, i!=j, represents the relative cost of retrieving 
     * entities of dimension i adjacent to entities of dimension j.
     * Diagonal entries, i==j, indicate whether or not handles to ALL
     * entities of that dimension are obtainable from calls that return
     * entity handles. This is always true by definition for i==j==0
     * as well as for i==j==2 in a 2D mesh and i==j==3 in a 3D mesh.
     * However, diagonal entries [1,1] for a 2D mesh and both [1,1]
     * and [2,2] for a 3D mesh indicate whether or not handles to ALL
     * intermediate dimensioned entities (ALL edges in a 2D mesh or
     * ALL edges and ALL faces in a 3D mesh) are obtainable from calls
     * that return entity handles. A value of iBase_AVAILABLE for a
     * diagonal entry indicates that handles are obtainable for ALL
     * entities of that dimension while a value of iBase_UNAVAILABLE
     * indicates that handles are not obtainable for ALL entities of that
     * dimension.
     * \param instance iMesh instance handle
     * \param *adjacency_table Pointer to array representing adjacency table
     *        returned from function
     * \param adjacency_table_allocated Pointer to allocated size of 
     *        adjacency table
     * \param adjacency_table_size Pointer to occupied size of 
     *        adjacency table
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAdjTable (iMesh_Instance instance,
                          /*out*/ int** adjacency_table,
                          /*inout*/ int* adjacency_table_allocated,
                          /*out*/ int* adjacency_table_size, 
                          /*out*/ int *err);

   /**\brief  Set the adjacency table as requested by the application
     *
     * Set the adjacency table as requested by the application. 
     * See iMesh_getAdjTable for a description of the meaning of the entries
     * in this table. This call to set the adjacency table allows an application
     * to request how it would like an implementation to behave when adjacent
     * and/or intermediate entities are queried. If the implementation is not
     * able to accommodate the requested query behavior associated with ANY
     * entry in the table, the call will fail and return an error of
     * iBase_NOT_SUPPORTED. Otherwise, the implementation is able to accommodate
     * the requested query behavior associated with ALL table entries and the
     * call will succeed. In either case, however, the implementation will
     * over-write the entries in the adjacency_table argument with the same
     * values that would be obtained in a succeeding call to iMesh_getAdjTable.
     * \param instance iMesh instance handle
     * \param adjacency_table Array representing adjacency table requested by application
     * \param adjacency_table_size Size of adj table array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setAdjTable (iMesh_Instance instance,
                          /*inout*/ int* adjacency_table,
                          /*in*/    int adjacency_table_size, 
                          /*out*/   int *err);

    /**\brief  Get the number of entities with the specified type in the instance or set
     *
     * Get the number of entities with the specified type in the instance 
     * or set.  If entity set handle is root set, return information for instance,
     * otherwise for set.  Value of entity type must be from the
     * iBase_EntityType enumeration.  If iBase_ALL_TYPES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entity requested
     * \param num_type Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfType(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_type,
                          /*out*/ int *num_type, 
                          /*out*/ int *err);

    /**\brief  Get the number of entities with the specified topology in the instance or set
     *
     * Get the number of entities with the specified topology in the instance 
     * or set.  If entity set handle is root set, return information for instance,
     * otherwise for set.  Value of entity topology must be from the
     * iMesh_EntityTopology enumeration.  If iMesh_ALL_TOPOLOGIES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_topology Topology of entity requested
     * \param num_topo Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfTopo(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_topology,
                          /*out*/ int *num_topo, 
                          /*out*/ int *err);

    /**\brief Return whether entity handles have changed since last reset or since
     *        instance construction
     *
     * Return whether entity handles have changed since last reset or since
     * instance construction.  If non-zero value is returned, it is not
     * guaranteed that a handle from before the last call to this function
     * represents the same entity as the same handle value does now.  If
     * doReset is non-zero, resets the starting point for this function.
     * \param instance iMesh instance handle
     * \param doReset Perform a reset on the starting point after which handles
     *        are invariant.
     * \param areHandlesInvariant Pointer to invariant flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_areEHValid(iMesh_Instance instance, 
                        /*in*/ int doReset,
                        /*out*/ int *areHandlesInvariant, 
                        /*out*/ int *err);


    /**\brief  Get entities of specific type and/or topology in set or instance
     *
     * Get entities of specific type and/or topology in set or instance.  All 
     * entities of a given type or topology are requested by specifying
     * iBase_ALL_TOPOLOGIES or iBase_ALL_TYPES, respectively.  Specified type
     * or topology must be a value in the iBase_EntityType or iBase_EntityTopology
     * enumeration, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entities being requested
     * \param entity_topology Topology of entities being requested
     * \param *entity_handles Pointer to array of entity handles returned 
     *        from function
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     *
     * Note: This function will fail and return an error if the caller passes a
     * combination of entity_type and entity_topology that are not consistent.
     */
  void iMesh_getEntities(iMesh_Instance instance,
                         /*in*/ const iBase_EntitySetHandle entity_set_handle,
                         /*in*/ const int entity_type,
                         /*in*/ const int entity_topology,
                         /*inout*/ iBase_EntityHandle** entity_handles,
                         /*inout*/ int* entity_handles_allocated,
                         /*out*/ int* entity_handles_size,
                         /*out*/ int *err);

    /**\brief  Get coordinates of specified vertices
     *
     * Get coordinates of specified vertices. Coordinates are returned
     * in the storage order indicated by the storage_order argument.
     * \param instance iMesh instance handle
     * \param vertex_handles Array of mesh vertex handles whose coordinates are
     *        being requested
     * \param vertex_handles_size Number of vertices in vertex_handles array
     * \param storage_order Requested storage order of returned coordinates
     * \param *coords Pointer to array of coordinates returned from function
     * \param *coords_allocated Pointer to allocated size of coords array
     * \param *coords_size Pointer to occupied size of coords array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getVtxArrCoords(iMesh_Instance instance,
                             /*in*/ const iBase_EntityHandle* vertex_handles,
                             /*in*/ const int vertex_handles_size,
                             /*in*/ const int storage_order,
                             /*inout*/ double** coords,
                             /*inout*/ int* coords_allocated,
                             /*out*/ int* coords_size, 
                             /*out*/ int *err);


    /**\brief Initialize an array iterator over specified entity type, topology, and 
     *        size
     *
     * Initialize an array iterator over specified entity type, topology, and 
     * size, for a specified set or instance.  Iterator returned can be used 
     * as input to functions returning entities for the iterator.  If all 
     * entities of a specified type and/or topology are to be iterated, 
     * specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  
     * Specified type or topology must be a value in the iBase_EntityType or 
     * iMesh_EntityTopology enumerations, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being iterated
     * \param requested_entity_type Type of entity to iterate
     * \param requested_entity_topology Topology of entity to iterate
     * \param requested_array_size Size of chunks of handles returned for each
     *        value of the iterator
     * \param entArr_iterator Pointer to iterator returned from function
     * \param *err Pointer to error type returned from function
     *
     * Note: This function will fail and return an error if the caller passes a
     * combination of entity_type and entity_topology that are not consistent.
     */
  void iMesh_initEntArrIter(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const int requested_entity_type,
                            /*in*/ const int requested_entity_topology,
                            /*in*/ const int requested_array_size,
                            /*out*/ iBase_EntityArrIterator* entArr_iterator,
                            /*out*/ int *err);

    /**\brief  Get entities contained in array iterator and increment iterator
     *
     * Get the entities corresponding to an array iterator (e.g. dereference
     * the array iterator), and increment the iterator. The dereferenced
     * value(s) are returned in entity_handles. If the iterator is at the
     * end of the iteration, the dereferenced value(s) are undefined and
     * has_data will be returned with a value of zero. Otherwise, has_data
     * will be returned with a non-zero value.
     * \param instance iMesh instance handle
     * \param entArr_iterator Iterator being queried
     * \param *entity_handles Pointer to array of entity handles contained in
     *        current value of iterator
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles 
     *        array
     * \param has_data Pointer to a flag indicating if the value(s) returned
              in entity_handles are valid. A non-zero value indicates the value(s)
              are valid. A zero value indicates the value(s) are NOT valid.
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNextEntArrIter(iMesh_Instance instance,
                               /*in*/ iBase_EntityArrIterator entArr_iterator,
                               /*inout*/ iBase_EntityHandle** entity_handles,
                               /*inout*/ int* entity_handles_allocated,
                               /*out*/ int* entity_handles_size,
                               /*out*/ int *has_data, 
                               /*out*/ int *err);


    /**\brief  Reset the array iterator
     *
     * Reset the array iterator
     * \param instance iMesh instance handle
     * \param entArr_iterator Iterator to reset
     * \param *err Pointer to error type returned from function
     */
  void iMesh_resetEntArrIter(iMesh_Instance instance,
                             /*in*/ iBase_EntityArrIterator entArr_iterator, 
                             /*out*/ int *err);


    /**\brief  Destroy the specified array iterator
     *
     * Destroy the specified array iterator
     * \param instance iMesh instance handle
     * \param entArr_iterator Iterator which gets destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_endEntArrIter(iMesh_Instance instance,
                           /*in*/ iBase_EntityArrIterator entArr_iterator, 
                           /*out*/ int *err);

    /**\brief  Get the entity topology for the specified entities
     *
     * Get the entity topology for the specified entities.  Topologies 
     * returned are values in the iMesh_EntityTopology enumeration.
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles being queried
     * \param entity_handles_size Number of entities in entity_handles array
     * \param *topology Pointer to array of entity topologies returned 
     *        from function
     * \param *topology_allocated Pointer to allocated size of topology array
     * \param *topology_size Pointer to occupied size of topology array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntArrTopo(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** topology,
                           /*inout*/ int* topology_allocated,
                           /*out*/ int* topology_size, 
                           /*out*/ int *err);


    /**\brief  Get the entity type for the specified entities
     *
     * Get the entity type for the specified entities.  Types
     * returned are values in the iBase_EntityType enumeration.
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles being queried
     * \param entity_handles_size Number of entities in entity_handles array
     * \param *type Pointer to array of types returned from function
     * \param *type_allocated Pointer to allocated size of type array
     * \param *type_size Pointer to occupied size of type array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntArrType(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** type,
                           /*inout*/ int* type_allocated,
                           /*out*/ int* type_size, 
                           /*out*/ int *err);


    /**\brief  Get entities of specified type adjacent to entities
     *
     * Get entities of specified type adjacent to entities.  Specified type
     * must be value in the iBase_EntityType enumeration.  \em offset(i) is
     * index of first entity in adjacentEntityHandles array adjacent to 
     * entity_handles[i].  More precisely, the entities adjacent to the
     * ith entity in entity_handles are found in adjacentEntityHandles
     * running from offset[i] to offset[i+1] - 1.  This implies that the
     * offset_size will be entity_handles_size + 1.
     *
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles being queried
     * \param entity_handles_size Number of entities in entity_handles array
     * \param entity_type_requested Type of adjacent entities requested
     * \param *adjacentEntityHandles Pointer to array of adjacentEntityHandles 
     *        returned from function. Note that the implicit INTERLEAVED storage
     *        order rule applies (see section ITAPS Storage Orders)
     * \param *adjacentEntityHandles_allocated Pointer to allocated size of 
     *        adjacentEntityHandles array
     * \param *adj_entity_handles_size Pointer to occupied size of 
     *        adjacentEntityHandles array
     * \param *offset Pointer to array of offsets returned from function
     * \param *offset_allocated Pointer to allocated size of offset array
     * \param *offset_size Pointer to occupied size of offset array
     * \param *err Pointer to error type returned from function
     *
     * Note 1: Because 'adjacent' as defined by the iMesh data model refers
     *         to those entities that bound another, the entities being queried
     *         here (in entity_handles arg) are NEVER ALSO returned in
     *         adjacentEntityHandles even if the entity_type_requested
     *         matches the entity type(s) in entity_handles.
     */
  void iMesh_getEntArrAdj(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const int entity_type_requested,
                          /*inout*/ iBase_EntityHandle** adjacentEntityHandles,
                          /*inout*/ int* adjacentEntityHandles_allocated,
                          /*out*/ int* adj_entity_handles_size,
                          /*inout*/ int** offset,
                          /*inout*/ int* offset_allocated,
                          /*out*/ int* offset_size,
                          /*out*/ int *err);

/**\brief Get "2nd order" adjacencies to an array of entities
 * Get "2nd order" adjacencies to an array of entities, that is, from each 
 * entity, through other entities of a specified "bridge" dimension, to other 
 * entities of another specified "to" dimension.
 *
 * Note 1:  If the "bridge" dimension is the same as the "to" dimension,
 *    the output will be empty (and an error code of
 *    iBase_INVALID_ARGUMENT returned).  If the type of a particular
 *    entity matches the "bridge" dimension, there will be no entities
 *    returned for that input entity.  This is consistent with the
 *    definition of adjacencies and the behavior of iMesh first
 *    adjacency calls. 
 * Note 2:  An entity will never be returned as a second adjacency of
 *    itself, on the grounds that this is the most likely expectation of
 *    applications, and that it is easier for an application to add the
 *    original entity to the returned data than to find and remove it.
 * Note 3: The entities adjacent to the ith entity in entity_handles are
 *    found in adj_entity_handles running from offset[i] to
 *    offset[i+1] - 1.  This implies that the offset_size will be
 *    entity_handles_size + 1. 
 *
 * \param instance iMesh instance for this call
 * \param entity_handles Entities from which adjacencies are requested
 * \param entity_handles_size Number of entities whose adjacencies are requested
 * \param bridge_entity_type  Type of bridge entity for 2nd order adjacencies
 * \param requested_entity_type Type of adjacent entities returned
 * \param adj_entity_handles Adjacent entities. Note that the implicit
 *        INTERLEAVED storage order rule applies
 *        (see section ITAPS Storage Orders)
 * \param adj_entity_handles_allocated Allocated size of returned array
 * \param adj_entity_handles_size Occupied size of returned array
 * \param offset Offset[i] is offset into adj_entity_handles of 2nd order 
 *        adjacencies of ith entity in entity_handles.  
 * \param offset_allocated Allocated size of offset array
 * \param offset_size Occupied size of offset array
 * \param err 
 */
  void iMesh_getEntArr2ndAdj( iMesh_Instance instance,
                              /*in*/ iBase_EntityHandle const* entity_handles,
                              /*in*/ int entity_handles_size,
                              /*in*/ int bridge_entity_type,
                              /*in*/ int requested_entity_type,
                              /*inout*/ iBase_EntityHandle** adj_entity_handles,
                              /*inout*/ int* adj_entity_handles_allocated,
                              /*out*/ int* adj_entity_handles_size,
                              /*inout*/ int** offset,
                              /*inout*/ int* offset_allocated,
                              /*out*/ int* offset_size,
                              /*out*/ int* err );

   /**\brief Get indexed representation of mesh or subset of mesh
    *
    * Given an entity set and optionally a type or topology, return:
    * - The entities in the set of the specified type or topology
    * - The entities adjacent to those entities with a specified
    *    type, as a list of unique handles.
    * - For each entity in the first list, the adjacent entities,
    *    specified as indices into the second list.
    *
    *\param entity_set_handle     The set of entities from which to query
    *\param entity_type_requester If not iBase_ALL_TYPES, act only on 
    *                             the subset of 'entity_set_handle' of the
    *                             specified type.
    *\param entity_topology_requester If not iMesh_ALL_TOPOLOGIES, act only
    *                             on the subset of 'entity_set_handle' with
    *                             the specified topology.
    *\param entity_type_requested The type of the adjacent entities to
    *                             return.
    *\param entity_handles        The handles of the (non-struct) subset of   
    *                             the entity set indicated by 
    *                             'entity_set_handle' and the optional type
    *                             and topology filtering arguments.
    *\param adj_entity_handles    The union of the unique entities of type 
    *                             'requested_entity_type' adjacent to each
    *                             entity in 'entity_handles'. Note that the
    *                             implicit INTERLEAVED storage order rule
    *                             applies (see section ITAPS Storage Orders)
    *\param adj_entity_indices    For each entity in 'entity_handles', the
    *                             adjacent entities of type
    *                             'entity_type_requested', specified as 
    *                             indices into 'adj_entity_handles'.  The
    *                             values are concatenated into a single   
    *                             array in the order of the entity handles 
    *                             in 'entity_handles'.
    *\param offset                For each entity in the corresponding 
    *                             position in 'entity_handles', the position
    *                             in 'adj_entity_indices' at which values
    *                             for that entity are stored.
    *
    * Note 1: Because 'adjacent' as defined by the iMesh data model refers
    *         to those entities that bound another, the entities being queried
    *         here (in entity_set_handle arg) are NEVER ALSO returned in
    *         adj_entity_handles even if the entity_type_requested
    *         matches the entity type(s) in entity_set_handle.
    * Note 2: The entities adjacent to the ith entity in entity_handles are
    *         found in adj_entity_handles running from offset[i] to
    *         offset[i+1] - 1.  This implies that the offset_size will
    *         be entity_handles_size + 1.  
    * Note 3: This function will fail and return an error if the caller passes
    *         a combination of entity_type and entity_topology that are
    *         not consistent. 
    */
  void iMesh_getAdjEntIndices(iMesh_Instance instance,
                      /*in*/    iBase_EntitySetHandle entity_set_handle,
                      /*in*/    int entity_type_requester,
                      /*in*/    int entity_topology_requester,
                      /*in*/    int entity_type_requested,
                      /*inout*/ iBase_EntityHandle** entity_handles,
                      /*inout*/ int* entity_handles_allocated,
                      /*out*/   int* entity_handles_size,
                      /*inout*/ iBase_EntityHandle** adj_entity_handles,
                      /*inout*/ int* adj_entity_handles_allocated,
                      /*out*/   int* adj_entity_handles_size,
                      /*inout*/ int** adj_entity_indices,
                      /*inout*/ int* adj_entity_indices_allocated,
                      /*out*/   int* adj_entity_indices_size,
                      /*inout*/ int** offset,
                      /*inout*/ int* offset_allocated,
                      /*out*/   int* offset_size,
                      /*out*/   int *err);

    /**\brief  Create an entity set
     *
     * Create an entity set, either ordered (isList=1) or unordered 
     * (isList=0).  Unordered entity sets can contain a given entity or 
     * set only once.
     * \param instance iMesh instance handle
     * \param isList If non-zero, an ordered list is created, otherwise an
     *        unordered set is created.
     * \param entity_set_created Entity set created by function
     * \param *err Pointer to error type returned from function
     */

  void iMesh_createEntSet(iMesh_Instance instance,
                          /*in*/ const int isList,
                          /*out*/ iBase_EntitySetHandle* entity_set_created,
                          /*out*/ int *err);

    /**\brief  Destroy an entity set
     *
     * Destroy an entity set
     * \param instance iMesh instance handle
     * \param entity_set Entity set to be destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_destroyEntSet(iMesh_Instance instance,
                           /*in*/ iBase_EntitySetHandle entity_set,
                           /*out*/ int *err);

    /**\brief  Return whether a specified set is ordered or unordered
     *
     * Return whether a specified set is ordered (*is_list=1) or 
     * unordered (*is_list=0)
     * \param instance iMesh instance handle
     * \param entity_set Entity set being queried
     * \param is_list Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isList(iMesh_Instance instance,
                    /*in*/ const iBase_EntitySetHandle entity_set,
                    /*out*/ int *is_list,
                    /*out*/ int *err);

    /**\brief  Get the number of entity sets contained in a set or interface
     *
     * Get the number of entity sets contained in a set or interface.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to contained set,
     *        not inclusive of the contained set. See notes on indirection
     *        in queries on set-inclusion structure.
     * \param num_sets Pointer to the number of sets returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumEntSets(iMesh_Instance instance,
                           /*in*/ const iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const int num_hops,
                           /*out*/ int *num_sets,
                           /*out*/ int *err);


    /**\brief  Get the entity sets contained in a set or interface
     *
     * Get the entity sets contained in a set or interface.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to contained set,
     *        not inclusive of the contained set. See notes on indirection
     *        in queries on set-inclusion structure.
     * \param *contained_set_handles Pointer to array of set handles returned
     *        from function
     * \param contained_set_handles_allocated Pointer to allocated length of
     *        contained_set_handles array
     * \param contained_set_handles_size Pointer to occupied length of
     *        contained_set_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSets(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set_handle,
                        /*in*/ const int num_hops,
                        /*inout*/ iBase_EntitySetHandle** contained_set_handles,
                        /*inout*/ int* contained_set_handles_allocated,
                        /*out*/ int* contained_set_handles_size,
                        /*out*/ int *err);

    /**\brief  Add an entity to a set
     *
     * Add an entity to a set
     * \param instance iMesh instance handle
     * \param entity_handle The entity being added
     * \param entity_set Pointer to the set being added to
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addEntToSet(iMesh_Instance instance,
                         /*in*/ iBase_EntityHandle entity_handle,
                         /*in*/ iBase_EntitySetHandle entity_set,
                         /*out*/ int *err);

    /**\brief  Remove an entity from a set
     *
     * Remove an entity from a set
     *
     * \param instance iMesh instance handle
     * \param entity_handle The entity being removed
     * \param entity_set Pointer to the set being removed from
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntFromSet(iMesh_Instance instance,
                           /*in*/ iBase_EntityHandle entity_handle,
                           /*in*/ iBase_EntitySetHandle entity_set,
                           /*out*/ int *err);


    /**\brief  Add an array of entities to a set
     *
     * Add an array of entities to a set
     * \param instance iMesh instance handle
     * \param entity_handles Array of entities being added
     * \param entity_handles_size Number of entities in entity_handles array
     * \param entity_set Pointer to the set being added to
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addEntArrToSet(iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ int entity_handles_size,
                            /*in*/ iBase_EntitySetHandle entity_set,
                            /*out*/ int *err);


    /**\brief  Remove an array of entities from a set
     *
     * Remove an array of entities from a set
     * \param instance iMesh instance handle
     * \param entity_handles Array of entities being remove
     * \param entity_handles_size Number of entities in entity_handles array
     * \param entity_set Pointer to the set being removed from
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntArrFromSet(iMesh_Instance instance,
                              /*in*/ const iBase_EntityHandle* entity_handles,
                              /*in*/ int entity_handles_size,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*out*/ int *err);


    /**\brief  Add an entity set to a set
     *
     * Add an entity set to a set
     * \param instance iMesh instance handle
     * \param entity_set_to_add The entity set being added
     * \param entity_set_handle Pointer to the set being added to
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addEntSet(iMesh_Instance instance,
                       /*in*/ iBase_EntitySetHandle entity_set_to_add,
                       /*in*/ iBase_EntitySetHandle entity_set_handle,
                       /*out*/ int *err);


    /**\brief  Remove an entity set from a set
     *
     * Remove an entity set from a set
     * \param instance iMesh instance handle
     * \param entity_set_to_remove The entity set being removed
     * \param entity_set_handle Pointer to the set being removed from
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntSet(iMesh_Instance instance,
                       /*in*/ iBase_EntitySetHandle entity_set_to_remove,
                       /*in*/ iBase_EntitySetHandle entity_set_handle,
                       /*out*/ int *err);

    /**\brief  Return whether an entity is contained in another set
     *
     * Return whether an entity is contained (*is_contained=1) or not 
     * contained (*is_contained=0) in another set
     * \param instance iMesh instance handle
     * \param containing_entity_set Entity set being queried
     * \param contained_entity Entity potentially contained in 
     *        containing_entity_set
     * \param is_contained Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isEntContained(iMesh_Instance instance,
                            /*in*/ iBase_EntitySetHandle containing_entity_set,
                            /*in*/ iBase_EntityHandle contained_entity,
                            /*out*/ int *is_contained,
                            /*out*/ int *err);

    /**\brief  Return whether entities are contained in a set
     *
     * Return whether each entity is contained in the set.
     * \param instance iMesh instance handle
     * \param containing_entity_set Entity set being queried
     * \param entity_handles List of entities for which to check containment.
     * \param is_contained One value for each input entity, 1 if contained
     *          in set, zero otherwise.
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isEntArrContained( iMesh_Instance instance,
                         /*in*/ iBase_EntitySetHandle containing_set,
                         /*in*/ const iBase_EntityHandle* entity_handles,
                         /*in*/ int num_entity_handles,
                      /*inout*/ int** is_contained,
                      /*inout*/ int* is_contained_allocated,
                        /*out*/ int* is_contained_size,
                        /*out*/ int* err );

    /**\brief  Return whether an entity set is contained in another set
     *
     * Return whether a set is contained (*is_contained=1) or not contained
     * (*is_contained=0) in another set
     * \param instance iMesh instance handle
     * \param containing_entity_set Entity set being queried
     * \param contained_entity_set Entity set potentially contained in 
     *        containing_entity_set
     * \param is_contained Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isEntSetContained(iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle containing_entity_set,
                               /*in*/ const iBase_EntitySetHandle contained_entity_set,
                               /*out*/ int *is_contained,
                               /*out*/ int *err);

    /**\brief  Add parent/child links between two sets
     *
     * Add parent/child links between two sets.  Makes parent point to child
     * and child point to parent.
     * \param instance iMesh instance handle
     * \param parent_entity_set Pointer to parent set
     * \param child_entity_set Pointer to child set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addPrntChld(iMesh_Instance instance,
                         /*in*/ iBase_EntitySetHandle parent_entity_set,
                         /*in*/ iBase_EntitySetHandle child_entity_set,
                         /*out*/ int *err);

    /**\brief  Remove parent/child links between two sets
     *
     * Remove parent/child links between two sets.
     * \param instance iMesh instance handle
     * \param parent_entity_set Pointer to parent set
     * \param child_entity_set Pointer to child set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvPrntChld(iMesh_Instance instance,
                         /*in*/ iBase_EntitySetHandle parent_entity_set,
                         /*in*/ iBase_EntitySetHandle child_entity_set,
                         /*out*/ int *err);

    /**\brief  Return whether two sets are related by parent/child links
     *
     * Return whether two sets are related (*is_child=1) or not (*is_child=0)
     * by parent/child links
     * \param instance iMesh instance handle
     * \param parent_entity_set Pointer to parent set
     * \param child_entity_set Pointer to child set
     * \param is_child Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isChildOf(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle parent_entity_set,
                       /*in*/ const iBase_EntitySetHandle child_entity_set,
                       /*out*/ int *is_child,
                       /*out*/ int *err);

    /**\brief  Get the number of child sets linked from a specified set
     *
     * Get the number of child sets linked from a specified set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to child set,
     *        not inclusive of the child set. See notes on indirection
     *        in parent-child structure.
     * \param num_child Pointer to number of children returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumChld(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        /*out*/ int *num_child,
                        /*out*/ int *err);

    /**\brief  Get the number of parent sets linked from a specified set
     *
     * Get the number of parent sets linked from a specified set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to parent set,
     *        not inclusive of the parent set. See notes on indirection
     *        in parent-child structure.
     * \param num_parent Pointer to number of parents returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumPrnt(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        /*out*/ int *num_parent,
                        /*out*/ int *err);

    /**\brief  Get the child sets linked from a specified set
     *
     * Get the child sets linked from a specified set.
     * \param instance iMesh instance handle
     * \param from_entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to child set,
     *        not inclusive of the child set. See notes on indirection
     *        in parent-child structure.
     * \param *entity_set_handles Pointer to array of child sets
     *        returned from function
     * \param *entity_set_handles_allocated Pointer to allocated size of 
     *        entity_set_handles array
     * \param *entity_set_handles_size Pointer to occupied size of 
     *        entity_set_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getChldn(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*inout*/ iBase_EntitySetHandle** entity_set_handles,
                      /*inout*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size,
                      /*out*/ int *err);

    /**\brief  Get the parent sets linked from a specified set
     *
     * Get the parent sets linked from a specified set.
     * \param instance iMesh instance handle
     * \param from_entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to parent set,
     *        not inclusive of the parent set. See notes on indirection
     *        in parent-child structure.
     * \param *entity_set_handles Pointer to array of parent sets
     *        returned from function
     * \param *entity_set_handles_allocated Pointer to allocated size of 
     *        entity_set_handles array
     * \param *entity_set_handles_size Pointer to occupied size of 
     *        entity_set_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getPrnts(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*inout*/ iBase_EntitySetHandle** entity_set_handles,
                      /*inout*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size,
                      /*out*/ int *err);

    /**\brief  Set coordinates for an array of vertices
     *
     * Set coordinates for an array of vertices.  Specified storage 
     * order must be either iBase_INTERLEAVED or iBase_BLOCKED, and 
     * indicates order of x, y, and z coordinates in coordinate array.
     * \param instance iMesh instance handle
     * \param vertex_handles Array of vertex handles
     * \param vertex_handles_size Number of vertex handles in array
     * \param storage_order Storage order of coordinates in coordinate array
     * \param new_coords Coordinate array
     * \param new_coords_size Size of coordinate array; should be 
     *        3*vertex_handles_size
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setVtxArrCoords(iMesh_Instance instance,
                             /*in*/ const iBase_EntityHandle* vertex_handles,
                             /*in*/ const int vertex_handles_size,
                             /*in*/ const int storage_order,
                             /*in*/ const double* new_coords,
                             /*in*/ const int new_coords_size,
                             /*out*/ int *err);


    /**\brief  Create an array of new vertices at specified coordinates
     *
     * Create an array of new vertices at specified coordinates.  Value of
     * storage_order must be either iBase_INTERLEAVED or iBase_BLOCKED.
     * \param instance iMesh instance handle
     * \param num_verts Number of new vertices to be created
     * \param storage_order Storage order of coordinates in new_coords array
     * \param new_coords Array of coordinates of new vertices
     * \param new_coords_size Number of coordinates in new_coords array, should
     *        be 3*num_verts
     * \param *new_vertex_handles Pointer to array of new vertex handles 
     *        returned from function
     * \param *new_vertex_handles_allocated Pointer to allocated size of 
     *        new_vertex_handles array
     * \param *new_vertex_handles_size Pointer to occupied size of 
     *        new_vertex_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createVtxArr(iMesh_Instance instance,
                          /*in*/ const int num_verts,
                          /*in*/ const int storage_order,
                          /*in*/ const double* new_coords,
                          /*in*/ const int new_coords_size,
                          /*inout*/ iBase_EntityHandle** new_vertex_handles,
                          /*inout*/ int* new_vertex_handles_allocated,
                          /*out*/ int* new_vertex_handles_size,
                          /*out*/ int *err);


    /**\brief  Create an array of new entities with specified lower-order topology
     *
     * Create an array of new entities with specified lower-order topology.  
     * Specified new_entity_topology must be value in iMesh_EntityTopology
     * enumeration.  Values return in status array must be values in the
     * iBase_CreationStatus enumeration.
     * \param instance iMesh instance handle
     * \param new_entity_topology Topology of created entity
     * \param lower_order_entity_handles Array of lower order entity handles
     *        used to construct new entities
     * \param lower_order_entity_handles_size Number of entities in array of 
     *        lower order entity handles
     * \param *new_entity_handles Pointer to array of new_entity_handles 
     *        returned from function
     * \param *new_entity_handles_allocated Pointer to allocated size of 
     *        new_entity_handles array
     * \param *new_entity_handles_size Pointer to occupied size of 
     *        new_entity_handles array
     * \param *status Pointer to array of creation status returned from 
     *        function
     * \param *status_allocated Pointer to allocated size of status array
     * \param *status_size Pointer to occupied size of status array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createEntArr(iMesh_Instance instance,
                          /*in*/ const int new_entity_topology,
                          /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                          /*in*/ const int lower_order_entity_handles_size,
                          /*inout*/ iBase_EntityHandle** new_entity_handles,
                          /*inout*/ int* new_entity_handles_allocated,
                          /*out*/ int* new_entity_handles_size,
                          /*inout*/ int** status,
                          /*inout*/ int* status_allocated,
                          /*out*/ int* status_size,
                          /*out*/ int *err);


    /**\brief  Delete specified entities
     *
     * Delete specified entities
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles to be deleted
     * \param entity_handles_size Number of entities in array to be deleted
     * \param *err Pointer to error type returned from function
     */
  void iMesh_deleteEntArr(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*out*/ int *err);


    /**\brief  Create a tag with specified name, size, and type
     *
     * Create a tag with specified name, size, and type.  Tag size is in
     * units of size of tag_type data types.  Value input for tag type must be 
     * value in iBase_TagType enumeration.
     * \param instance iMesh instance handle
     * \param tag_name Character string indicating tag name
     * \param tag_size Size of each tag value, in units of number of tag_type 
     *        entities
     * \param tag_type Data type for data stored in this tag
     * \param tag_handle Pointer to tag handle returned from function
     * \param *err Pointer to error type returned from function
     * \param tag_name_len Length of tag name string
     */
  void iMesh_createTag(iMesh_Instance instance,
                       /*in*/ const char* tag_name,
                       /*in*/ const int tag_size,
                       /*in*/ const int tag_type,
                       /*out*/ iBase_TagHandle* tag_handle, 
                       /*out*/ int *err,
                       /*in*/ const int tag_name_len);


    /**\brief  Destroy a tag
     *
     * Destroy a tag.  If forced is non-zero and entities still have values
     * set for this tag, tag is deleted anyway and those values disappear,
     * otherwise tag is not deleted.
     * \param instance iMesh instance handle
     * \param tag_handle Handle of tag to be deleted
     * \param forced If non-zero, delete the tag even if entities have values
     *        set for that tag
     * \param *err Pointer to error type returned from function
     */
  void iMesh_destroyTag(iMesh_Instance instance,
                        /*in*/ iBase_TagHandle tag_handle,
                        /*in*/ const int forced,
                        /*out*/ int *err);

    /**\brief  Get the name for a given tag handle
     *
     * Get the name for a given tag handle
     * \param instance iMesh instance handle
     * \param tag_handle Tag handle being queried
     * \param name Pointer to character string to store name returned from 
     *        function
     * \param *err Pointer to error type returned from function
     * \param name_len Length of character string input to function
     */
  void iMesh_getTagName(iMesh_Instance instance,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*inout*/ char *name,
                        /*out*/ int *err,
                        /*in*/ int name_len);

    /**\brief  Get size of a tag in units of numbers of tag data type
     *
     * Get size of a tag in units of numbers of tag data type
     * \param instance iMesh instance handle
     * \param tag_handle Handle of tag being queried
     * \param tag_size Pointer to tag size returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getTagSizeValues(iMesh_Instance instance,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ int *tag_size,
                              /*out*/ int *err);

    /**\brief  Get size of a tag in units of bytes
     *
     * Get size of a tag in units of bytes
     * \param instance iMesh instance handle
     * \param tag_handle Handle of tag being queried
     * \param tag_size Pointer to tag size returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getTagSizeBytes(iMesh_Instance instance,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*out*/ int *tag_size,
                             /*out*/ int *err);

    /**\brief  Get a the handle of an existing tag with the specified name
     *
     * Get a the handle of an existing tag with the specified name
     * \param instance iMesh instance handle
     * \param tag_name Name of tag being queried
     * \param tag_handle Pointer to tag handle returned from function
     * \param *err Pointer to error type returned from function
     * \param tag_name_len Length of tag name string
     */
  void iMesh_getTagHandle(iMesh_Instance instance,
                          /*in*/ const char* tag_name,
                          /*out*/ iBase_TagHandle *tag_handle, 
                          /*out*/ int *err,
                          int tag_name_len);

    /**\brief  Get the data type of the specified tag handle
     *
     * Get the data type of the specified tag handle.  Tag type is a value in
     * the iBase_TagType enumeration.
     * \param instance iMesh instance handle
     * \param tag_handle Handle for the tag being queried
     * \param tag_type Pointer to tag type returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getTagType(iMesh_Instance instance,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ int *tag_type,
                        /*out*/ int *err);

    /**\brief  Set a tag value of arbitrary type on an entity set
     *
     * Set a tag value of arbitrary type on an entity set. The tag data
     * is passed as void*. tag_value_size specifies the size of the memory
     * pointed to by tag_value in terms of bytes. Applications are free to
     * use this function to set data of any type, not just iBase_BYTES.
     * However, in all cases, the size specified by tag_value_size is
     * always in terms of bytes.
     *
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Pointer to tag data being set on entity set
     * \param tag_value_size Size in bytes of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetData(iMesh_Instance instance,
                           /*in*/ iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const void* tag_value,
                           /*in*/ const int tag_value_size,
                           /*out*/ int *err);


    /**\brief  Set a tag value of integer type on an entity set
     *
     * Set a tag value of integer type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Tag value being set on entity set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetIntData(iMesh_Instance instance,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*in*/ const int tag_value,
                              /*out*/ int *err);


    /**\brief  Set a tag value of double type on an entity set
     *
     * Set a tag value of double type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Tag value being set on entity set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetDblData(iMesh_Instance instance,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*in*/ const double tag_value,
                              /*out*/ int *err);


    /**\brief  Set a tag value of entity handle type on an entity set
     *
     * Set a tag value of entity handle type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Tag value being set on entity set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetEHData(iMesh_Instance instance,
                             /*in*/ iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*in*/ const iBase_EntityHandle tag_value,
                             /*out*/ int *err);

    /**\brief  Set a tag value of entity set handle type on an
     *         entity set
     *
     * Set a tag value of entity set handle type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Tag value being set on entity set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetESHData(iMesh_Instance instance,
                             /*in*/ iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*in*/ const iBase_EntitySetHandle tag_value,
                             /*out*/ int *err);


    /**\brief  Get the value of a tag of arbitrary type on an entity set
     *
     * Get the value of a tag of arbitrary type on an entity set.  Tag data 
     * is returned back as void*. tag_value_size specifies the size of the
     * memory pointed to by tag_value in terms of bytes. Applications may
     * use this function to get data of any type, not just iBase_BYTES.
     * However because this function supports data of arbitrary type,
     * in all cases the size specified by tag_value_size is always in terms
     * of bytes.
     *
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *tag_value Pointer to tag data array being queried
     * \param *tag_value_allocated Pointer to tag data array allocated size
     * \param *tag_value_size Pointer to occupied size in bytes of tag data
     *        array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetData(iMesh_Instance instance,
                           /*in*/ const iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ void** tag_value,
                           /*inout*/ int* tag_value_allocated,
                           /*out*/ int* tag_value_size,
                           /*out*/ int *err);

    /**\brief  Get the value of a tag of integer type on an entity set
     *
     * Get the value of a tag of integer type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetIntData(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ int *out_data,
                              /*out*/ int *err);

    /**\brief  Get the value of a tag of double type on an entity set
     *
     * Get the value of a tag of double type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetDblData(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ double *out_data,
                              /*out*/ int *err);

    /**\brief  Get the value of a tag of entity handle type on an entity set
     *
     * Get the value of a tag of entity handle type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetEHData(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*out*/ iBase_EntityHandle *out_data,
                             /*out*/ int *err);

    /**\brief  Get the value of a tag of entity set handle type on an
     *         entity set
     *
     * Get the value of a tag of entity set handle type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetESHData(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*out*/ iBase_EntitySetHandle *out_data,
                             /*out*/ int *err);

    /**\brief  Get all the tags associated with a specified entity set
     *
     * Get all the tags associated with a specified entity set
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity being queried
     * \param *tag_handles Pointer to array of tag_handles returned from 
     *        function
     * \param *tag_handles_allocated Pointer to allocated size of tag_handles 
     *        array
     * \param *tag_handles_size Pointer to occupied size of tag_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAllEntSetTags(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set_handle,
                              /*inout*/ iBase_TagHandle** tag_handles,
                              /*inout*/ int* tag_handles_allocated,
                              /*out*/ int* tag_handles_size,
                              /*out*/ int *err);

    /**\brief  Remove a tag value from an entity set
     *
     * Remove a tag value from an entity set
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set from which tag is being removed
     * \param tag_handle Tag handle of tag being removed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntSetTag(iMesh_Instance instance,
                          /*in*/ iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*out*/ int *err);

    /**\brief  Set coordinates for a vertex
     *
     * Set coordinates for a vertex.
     * \param instance iMesh instance handle
     * \param vertex_handle vertex handle being set
     * \param x x coordinate being set
     * \param y y coordinate being set
     * \param z z coordinate being set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setVtxCoord(iMesh_Instance instance,
                         /*in*/ iBase_EntityHandle vertex_handle,
                         /*in*/ const double x, /*in*/ const double y,
                         /*in*/ const double z,
                         /*out*/ int *err);

    /**\brief  Create a new vertex at specified coordinates
     *
     * Create a new vertex at specified coordinates.
     * \param instance iMesh instance handle
     * \param x x coordinate of new vertex
     * \param y y coordinate of new vertex
     * \param z z coordinate of new vertex
     * \param new_vertex_handle Pointer to new vertex handles returned from 
     *        function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createVtx(iMesh_Instance instance,
                       /*in*/ const double x, /*in*/ const double y,
                       /*in*/ const double z,
                       /*out*/ iBase_EntityHandle* new_vertex_handle,
                       /*out*/ int *err);

    /**\brief  Create a new entity with specified lower-order topology
     *
     * Create a new entity with specified lower-order topology.  
     * Specified new_entity_topology must be value in iMesh_EntityTopology
     * enumeration.  Value returned as status must be a value in the
     * iBase_CreationStatus enumeration.
     * \param instance iMesh instance handle
     * \param new_entity_topology Topology of created entity
     * \param lower_order_entity_handles Array of lower order entity handles
     *        used to construct new entity
     * \param lower_order_entity_handles_size Number of entities in array of 
     *        lower order entity handles
     * \param new_entity_handle Pointer to new entity handle returned from 
     *        function
     * \param status Pointer to creation status returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createEnt(iMesh_Instance instance,
                       /*in*/ const int new_entity_topology,
                       /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                       /*in*/ const int lower_order_entity_handles_size,
                       /*out*/ iBase_EntityHandle* new_entity_handle,
                       /*out*/ int* status,
                       /*out*/ int *err);

    /**\brief  Delete specified entity
     *
     * Delete specified entity
     * \param instance iMesh instance handle
     * \param entity_handle Entity to be deleted
     * \param *err Pointer to error type returned from function
     */
  void iMesh_deleteEnt(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       /*out*/ int *err);

    /**\brief  Get tag values of arbitrary type for an array of entities
     *
     * Get tag values of arbitrary type for an array of entities.  Tag data 
     * is returned as void*. tag_values_size specifies the size of the
     * memory pointed to by tag_values in terms of bytes. Applications may
     * use this function to get data of any type, not just iBase_BYTES.
     * However, because this function supports data of arbitrary type, in
     * all cases the size specified by tag_values_size always in terms of
     * bytes.
     *
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_values Pointer to tag data array being returned from 
     *        function. Note that the implicit INTERLEAVED storage
     *        order rule applies (see section ITAPS Storage Orders)
     * \param tag_values_allocated Pointer to allocated size of tag data array
     * \param tag_values_size Pointer to occupied size in bytes of tag data
     *        array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getArrData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle* entity_handles,
                        /*in*/ const int entity_handles_size,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*inout*/ void** tag_values,
                        /*inout*/int* tag_values_allocated,
                        /*out*/ int* tag_values_size,
                        /*out*/ int *err);

    /**\brief  Get tag values of integer type for an array of entities
     *
     * Get tag values of integer type for an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_values Pointer to tag data array being returned from 
     *        function. Note that the implicit INTERLEAVED storage
     *        order rule applies (see section ITAPS Storage Orders)
     * \param tag_values_allocated Pointer to allocated size of tag data array
     * \param tag_values_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getIntArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ int** tag_values,
                           /*inout*/ int* tag_values_allocated,
                           /*out*/ int* tag_values_size,
                           /*out*/ int *err);

    /**\brief  Get tag values of double type for an array of entities
     *
     * Get tag values of double type for an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_values Pointer to tag data array being returned from 
     *        function. Note that the implicit INTERLEAVED storage
     *        order rule applies (see section ITAPS Storage Orders)
     * \param tag_values_allocated Pointer to allocated size of tag data array
     * \param tag_values_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getDblArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ double** tag_values,
                           /*inout*/ int* tag_values_allocated,
                           /*out*/ int* tag_values_size,
                           /*out*/ int *err);

    /**\brief  Get tag values of entity handle type for an array of entities
     *
     * Get tag values of entity handle type for an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_value Pointer to tag data array being returned from 
     *        function. Note that the implicit INTERLEAVED storage
     *        order rule applies (see section ITAPS Storage Orders)
     * \param tag_value_allocated Pointer to allocated size of tag data array
     * \param tag_value_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEHArrData(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*inout*/ iBase_EntityHandle** tag_value,
                          /*inout*/ int* tag_value_allocated,
                          /*out*/ int* tag_value_size,
                          /*out*/ int *err);

    /**\brief  Get tag values of entity set handle type for an array of
     *         entities
     *
     * Get tag values of entity set handle type for an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_value Pointer to tag data array being returned from 
     *        function. Note that the implicit INTERLEAVED storage
     *        order rule applies (see section ITAPS Storage Orders)
     * \param tag_value_allocated Pointer to allocated size of tag data array
     * \param tag_value_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getESHArrData(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*inout*/ iBase_EntitySetHandle** tag_value,
                          /*inout*/ int* tag_value_allocated,
                          /*out*/ int* tag_value_size,
                          /*out*/ int *err);

    /**\brief  Set tag values of arbitrary type on an array of entities
     *
     * Set tag values of arbitrary type on an array of entities.  Tag data
     * is passed as void*. tag_values_size specifies the size of the
     * memory pointed to by tag_values in terms of bytes. Applications may
     * use this function to set data of any type, not just iBase_BYTES.
     * However, because this function supports data of arbitrary type, in all
     * cases the size specified by tag_values_size is always in terms of
     * bytes.
     *
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entity. Note that
     *        the implicit INTERLEAVED storage order rule applies (see section
     *        ITAPS Storage Orders)
     * \param tag_values_size Size in bytes of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setArrData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle* entity_handles,
                        /*in*/ const int entity_handles_size,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const void* tag_values,
                        /*in*/ const int tag_values_size,
                        /*out*/ int *err);

    /**\brief  Set tag values of integer type on an array of entities
     *
     * Set tag values of integer type on an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entities. Note
     *        that the implicit INTERLEAVED storage order rule applies
     *        (see section ITAPS Storage Orders)
     * \param tag_values_size Size in total number of integers of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setIntArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const int* tag_values,
                           /*in*/ const int tag_values_size,
                           /*out*/ int *err);

    /**\brief  Set tag values of double type on an array of entities
     *
     * Set tag values of double type on an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entities. Note
     *        that the implicit INTERLEAVED storage order rule applies
     *        (see section ITAPS Storage Orders)
     * \param tag_values_size Size in total number of doubles of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setDblArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const double* tag_values,
                           /*in*/ const int tag_values_size,
                           /*out*/ int *err);

    /**\brief  Set tag values of entity handle type on an array of entities
     *
     * Set tag values of entity handle type on an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entities. Note
     *        that the implicit INTERLEAVED storage order rule applies
     *        (see section ITAPS Storage Orders)
     * \param tag_values_size Size in total number of entity handles of tag 
     *        data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEHArrData(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*in*/ const iBase_EntityHandle* tag_values,
                          /*in*/ const int tag_values_size,
                          /*out*/ int *err);

    /**\brief  Set tag values of entity set handle type on an array of
     *         entities
     *
     * Set tag values of entity set handle type on an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entities. Note
     *        that the implicit INTERLEAVED storage order rule applies
     *        (see section ITAPS Storage Orders)
     * \param tag_values_size Size in total number of entity handles of tag 
     *        data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setESHArrData(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*in*/ const iBase_EntitySetHandle* tag_values,
                          /*in*/ const int tag_values_size,
                          /*out*/ int *err);

    /**\brief  Remove a tag value from an array of entities
     *
     * Remove a tag value from an array of entities
     * \param instance iMesh instance handle
     * \param entity_handles Entity from which tag is being removed
     * \param entity_handles_size Number of entities in entity array
     * \param tag_handle Tag handle of tag being removed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvArrTag(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle* entity_handles,
                       /*in*/ const int entity_handles_size,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*out*/ int *err);

    /**\brief  Get the value of a tag of arbitrary type on an entity
     *
     * Get the value of a tag of arbitrary type on an entity.  Tag data 
     * is passed back as void*. tag_value_size specifies the size of the
     * memory pointed to by tag_value in terms of bytes. Applications may
     * use this function to get data of any type, not just iBase_BYTES.
     * However, because this function supports arbitrary type, in all
     * cases the size specified by tag_value_size is always in terms of
     * bytes.
     *
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *tag_value Pointer to tag data array being queried
     * \param *tag_value_allocated Pointer to tag data array allocated size
     * \param *tag_value_size Pointer to occupied size in bytes of tag data
     *        array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getData(iMesh_Instance instance,
                     /*in*/ const iBase_EntityHandle entity_handle,
                     /*in*/ const iBase_TagHandle tag_handle,
                     /*inout*/ void** tag_value,
                     /*inout*/ int *tag_value_allocated,
                     /*out*/ int *tag_value_size,
                     /*out*/ int *err);

    /**\brief  Get the value of a tag of integer type on an entity
     *
     * Get the value of a tag of integer type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getIntData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ int *out_data,
                        /*out*/ int *err);

    /**\brief  Get the value of a tag of double type on an entity
     *
     * Get the value of a tag of double type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getDblData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ double *out_data,
                        /*out*/ int *err);

    /**\brief  Get the value of a tag of entity handle type on an entity
     *
     * Get the value of a tag of entity handle type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEHData(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*out*/ iBase_EntityHandle *out_data,
                       /*out*/ int *err);

    /**\brief  Get the value of a tag of entity set handle type on an
     *         entity
     *
     * Get the value of a tag of entity set handle type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getESHData(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*out*/ iBase_EntitySetHandle *out_data,
                       /*out*/ int *err);

    /**\brief  Set a tag value of arbitrary type on an entity
     *
     * Set a tag value of arbitrary type on an entity.  Tag data
     * is passed as void*. tag_value_size specifies the size of the
     * memory pointed to by tag_value in terms of bytes. Applications may
     * use this function to set data of any type, not just iBase_BYTES.
     * However, because this function supports data of arbitrary type, in
     * all cases the size specified by tag_value_size is always in terms
     * of bytes.
     *
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Pointer to tag data being set on entity
     * \param tag_value_size Size in bytes of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setData(iMesh_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ const iBase_TagHandle tag_handle,
                     /*in*/ const void* tag_value,
                     /*in*/ const int tag_value_size,
                     /*out*/ int *err);

    /**\brief  Set a tag value of integer type on an entity
     *
     * Set a tag value of integer type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Tag value being set on entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setIntData(iMesh_Instance instance,
                        /*in*/ iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const int tag_value,
                        /*out*/ int *err);

    /**\brief  Set a tag value of double type on an entity
     *
     * Set a tag value of double type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Tag value being set on entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setDblData(iMesh_Instance instance,

                        /*in*/ iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const double tag_value,
                        /*out*/ int *err);

    /**\brief  Set a tag value of entity handle type on an entity
     *
     * Set a tag value of entity handle type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Tag value being set on entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEHData(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*in*/ const iBase_EntityHandle tag_value,
                       /*out*/ int *err);

    /**\brief  Set a tag value of entity set handle type on an entity
     *
     * Set a tag value of entity set handle type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Tag value being set on entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setESHData(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*in*/ const iBase_EntitySetHandle tag_value,
                       /*out*/ int *err);

    /**\brief  Get all the tags associated with a specified entity handle
     *
     * Get all the tags associated with a specified entity handle
     * \param instance iMesh instance handle
     * \param entity_handle Entity being queried
     * \param *tag_handles Pointer to array of tag_handles returned from 
     *        function
     * \param *tag_handles_allocated Pointer to allocated size of tag_handles 
     *        array
     * \param *tag_handles_size Pointer to occupied size of tag_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAllTags(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*inout*/ iBase_TagHandle** tag_handles,
                        /*inout*/ int* tag_handles_allocated,
                        /*out*/ int* tag_handles_size,
                        /*out*/ int *err);

    /**\brief  Remove a tag value from an entity
     *
     * Remove a tag value from an entity
     * \param instance iMesh instance handle
     * \param entity_handle Entity from which tag is being removed
     * \param tag_handle Tag handle of tag being removed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvTag(iMesh_Instance instance,
                    /*in*/ iBase_EntityHandle entity_handle,
                    /*in*/ const iBase_TagHandle tag_handle,
                    /*out*/ int *err);

    /**\brief  Initialize an iterator over specified entity type, topology, and size
     *
     * Initialize an iterator over specified entity type, topology, and size,
     * for a specified set or instance.  Iterator returned can be used as input
     * to functions returning the entity for the iterator.  If all entities of 
     * a specified type and/or topology are to be iterated, specify 
     * iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  Specified type 
     * or topology must be a value in the iBase_EntityType or 
     * iMesh_EntityTopology enumerations, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being iterated
     * \param requested_entity_type Type of entity to iterate
     * \param requested_entity_topology Topology of entity to iterate
     * \param entity_iterator Pointer to iterator returned from function
     * \param *err Pointer to error type returned from function
     *
     * Note: This function will fail and return an error if the caller passes a
     * combination of entity_type and entity_topology that are not consistent.
     */
  void iMesh_initEntIter(iMesh_Instance instance,
                         /*in*/ const iBase_EntitySetHandle entity_set_handle,
                         /*in*/ const int requested_entity_type,
                         /*in*/ const int requested_entity_topology,
                         /*out*/ iBase_EntityIterator* entity_iterator,
                         /*out*/ int *err);

    /**\brief  Get entity corresponding to an iterator and increment iterator
     *
     * Get the entity corresponding to an iterator (that is, dereference
     * the iterator), and increment the iterator. The dereferenced value
     * is returned in 'entity_handle'. If the iterator is at the end of the
     * iteration, the dereferenced value will be undefined and 'has_data'
     * will have a value of zero. Otherwise, 'has_data' will have a non-zero
     * value.
     * \param instance iMesh instance handle
     * \param entity_iterator Iterator being queried
     * \param entity_handle Pointer to an entity handle corresponding to the
     *        current value of iterator just prior to the call.
     * \param has_data Pointer to a flag indicating if the value returned in
              entity_handle is valid. A non-zero value indicates the value is
              valid. A zero value indicates the value is NOT valid.
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNextEntIter(iMesh_Instance instance,
                            /*in*/ iBase_EntityIterator entity_iterator,
                            /*out*/ iBase_EntityHandle* entity_handle,
                            /*out*/ int *has_data,
                            /*out*/ int *err);

    /**\brief  Reset the iterator
     *
     * Reset the iterator
     * \param instance iMesh instance handle
     * \param entity_iterator Iterator to reset
     * \param *err Pointer to error type returned from function
     */
  void iMesh_resetEntIter(iMesh_Instance instance,
                          /*in*/ iBase_EntityIterator entity_iterator,
                          /*out*/ int *err);

    /**\brief  Destroy the specified iterator
     *
     * Destroy the specified iterator
     * \param instance iMesh instance handle
     * \param entity_iterator Iterator which gets destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_endEntIter(iMesh_Instance instance,
                        /*in*/ iBase_EntityIterator entity_iterator,
                        /*out*/ int *err);

    /**\brief  Get the entity topology for the specified entity
     *
     * Get the entity topology for the specified entity.  Topology
     * returned is a value in the iMesh_EntityTopology enumeration.
     * \param instance iMesh instance handle
     * \param entity_handle Entity handle being queried
     * \param *out_topo Pointer to entity topology returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntTopo(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*out*/ int *out_topo,
                        /*out*/ int *err);

    /**\brief  Get the entity type for the specified entity
     *
     * Get the entity type for the specified entity.  Type returned is a value
     * in the iBase_EntityType enumeration.
     * \param instance iMesh instance handle
     * \param entity_handle Entity handle being queried
     * \param *out_type Pointer to entity type returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntType(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*out*/ int *out_type,
                        /*out*/ int *err);

    /**\brief  Get coordinates of specified vertex
     *
     * Get coordinates of specified vertex.
     * \param instance iMesh instance handle
     * \param vertex_handle Mesh vertex being queried
     * \param *x Pointer to x coordinate returned from function
     * \param *y Pointer to y coordinate returned from function
     * \param *z Pointer to z coordinate returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getVtxCoord(iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle vertex_handle,
                         /*out*/ double *x, /*out*/ double *y, /*out*/ double *z,
                         /*out*/ int *err);

    /**\brief  Get entities of specified type adjacent to an entity
     *
     * Get entities of specified type adjacent to an entity.  Specified type
     * must be value in the iBase_EntityType enumeration.
     * \param instance iMesh instance handle
     * \param entity_handle Entity handle being queried
     * \param entity_type_requested Type of adjacent entities requested
     * \param *adj_entity_handles Pointer to array of adjacent entities
     *        returned from function
     * \param *adj_entity_handles_allocated Pointer to allocated size of 
     *        adj_entity_handles array
     * \param *adj_entity_handles_size Pointer to occupied size of 
     *        adj_entity_handles array
     * \param *err Pointer to error type returned from function
     *
     * Note 1: Because 'adjacent' as defined by the iMesh data model refers
     *         to those entities that bound another, the entity being queried
     *         here (in entity_handle arg) is NEVER ALSO returned in
     *         adj_entity_handles even if the entity_type_requested
     *         matches the entity type in entity_handle.
     */
  void iMesh_getEntAdj(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const int entity_type_requested,
                       /*inout*/ iBase_EntityHandle** adj_entity_handles,
                       /*inout*/ int* adj_entity_handles_allocated,
                       /*out*/ int* adj_entity_handles_size,
                       /*out*/ int *err);

/**\brief Get "2nd order" adjacencies to an entity
 * Get "2nd order" adjacencies to an entity, that is, from an entity, through
 * other entities of a specified "bridge" dimension, to other entities of another 
 * specified "to" dimension.
 *
 * Note 1: If the "bridge" dimension is the same as the "to" dimension
 *    or the dimension of the input entity, the output will be empty
 *    (and an error code of iBase_INVALID_ARGUMENT returned).  This is
 *    consistent with the definition of adjacencies and the behavior of
 *    iMesh first adjacency calls.
 * Note 2: An entity will never be returned as a second adjacency of
 *    itself, on the grounds that this is the most likely expectation of
 *    applications, and that it is easier for an application to add the
 *    original entity to the returned data than to find and remove it.
 *
 * \param instance iMesh instance for this call
 * \param entity_handle Entity from which adjacencies are requested
 * \param bridge_entity_type  Type of bridge entity for 2nd order adjacencies
 * \param requested_entity_type Type of adjacent entities returned
 * \param adjacent_entities Adjacent entities
 * \param adjacent_entities_allocated Allocated size of returned array
 * \param adjacent_entities_size Occupied size of returned array
 * \param err 
 */
  void iMesh_getEnt2ndAdj( iMesh_Instance instance,
                           /*in*/ iBase_EntityHandle entity_handle,
                           /*in*/ int bridge_entity_type,
                           /*in*/ int requested_entity_type,
                           /*inout*/ iBase_EntityHandle** adjacent_entities,
                           /*inout*/ int* adjacent_entities_allocated,
                           /*out*/ int* adjacent_entities_size,
                           /*out*/ int* err );

    /**\brief  Subtract contents of one entity set from another
     *
     * Subtract contents of one entity set from another
     * \param instance iMesh instance handle
     * \param entity_set_1 Entity set from which other set is being subtracted
     * \param entity_set_2 Entity set being subtracted from other set
     * \param result_entity_set Pointer to entity set returned from function
     * \param *err Pointer to error type returned from function
     */

  void iMesh_subtract(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle entity_set_1,
                      /*in*/ const iBase_EntitySetHandle entity_set_2,
                      /*out*/ iBase_EntitySetHandle* result_entity_set,
                      /*out*/ int *err);

    /**\brief  Intersect contents of one entity set with another
     *
     * Intersect contents of one entity set with another
     * \param instance iMesh instance handle
     * \param entity_set_1 Entity set being intersected with another
     * \param entity_set_2 Entity set being intersected with another
     * \param result_entity_set Pointer to entity set returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_intersect(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle entity_set_1,
                       /*in*/ const iBase_EntitySetHandle entity_set_2,
                       /*out*/ iBase_EntitySetHandle* result_entity_set,
                       /*out*/ int *err);

    /**\brief  Unite contents of one entity set with another
     *
     * Unite contents of one entity set with another
     * \param instance iMesh instance handle
     * \param entity_set_1 Entity set being united with another
     * \param entity_set_2 Entity set being united with another
     * \param result_entity_set Pointer to entity set returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_unite(iMesh_Instance instance,
                   /*in*/ const iBase_EntitySetHandle entity_set_1,
                   /*in*/ const iBase_EntitySetHandle entity_set_2,
                   /*out*/ iBase_EntitySetHandle* result_entity_set,
                   /*out*/ int *err);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ifndef _ITAPS_iMesh */
