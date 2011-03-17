#ifndef _ITAPS_iBase
#define _ITAPS_iBase

#ifdef __cplusplus

extern "C"
{
#endif

#define IBASE_MINENUM(enumName) enumName ## _MIN
#define IBASE_MAXENUM(enumName) enumName ## _MAX
#define IBASE_NUMENUM(enumName) ((int)IBASE_MAXENUM(enumName) - (int)IBASE_MINENUM(enumName) + 1)
#define IBASE_INCENUM(enumName,I) (I = (enum enumName)((int)I+1))

    /*==========================================================
     * TYPEDEF'S
     *==========================================================
     */
  typedef void* iBase_Instance;
  typedef struct iBase_EntityHandle_Private* iBase_EntityHandle;
  typedef struct iBase_EntitySetHandle_Private* iBase_EntitySetHandle;
  typedef struct iBase_TagHandle_Private* iBase_TagHandle;
  typedef struct iBase_EntityIterator_Private* iBase_EntityIterator;
  typedef struct iBase_EntityArrIterator_Private* iBase_EntityArrIterator;


    /*==========================================================
     * ENTITYTYPE ENUMERATION
     *==========================================================
     */
  enum iBase_EntityType {
    iBase_EntityType_MIN = 0,
    iBase_VERTEX = iBase_EntityType_MIN,
    iBase_EDGE,
    iBase_FACE,
    iBase_REGION,
    iBase_ALL_TYPES,
    iBase_EntityType_MAX = iBase_ALL_TYPES
  };

    /*==========================================================
     * ADJACENCYCOST ENUMERATION
     *==========================================================
     */
  enum iBase_AdjacencyCost {
    iBase_AdjacencyCost_MIN = 0,
    iBase_UNAVAILABLE = iBase_AdjacencyCost_MIN, /**< Adjacency information not supported */
    iBase_ALL_ORDER_1,              /**< No more than local mesh traversal required (i!=j) */
    iBase_ALL_ORDER_LOGN,           /**< Global tree search (i!=j) */
    iBase_ALL_ORDER_N,              /**< Global exhaustive search (i!=j) */
    iBase_SOME_ORDER_1,             /**< Only some adjacency info, local (i!=j) */
    iBase_SOME_ORDER_LOGN,          /**< Only some adjacency info, tree (i!=j) */
    iBase_SOME_ORDER_N,             /**< Only some adjacency info, exhaustive (i!=j) */
    iBase_AVAILABLE,                /**< ALL (intermediate) entities available. (i==j) */
    iBase_AdjacencyCost_MAX = iBase_AVAILABLE
  };

    /*==========================================================
     * CREATIONSTATUS ENUMERATION
     *==========================================================
     */
  enum iBase_CreationStatus {
    iBase_CreationStatus_MIN = 0,
    iBase_NEW = iBase_CreationStatus_MIN,
    iBase_ALREADY_EXISTED,
    iBase_CREATED_DUPLICATE,
    iBase_CREATION_FAILED,
    iBase_CreationStatus_MAX = iBase_CREATION_FAILED
  };

    /*==========================================================
     * ERRORACTIONS ENUMERATION
     *==========================================================
     */
  enum iBase_ErrorActions {
    iBase_ErrorActions_MIN = 0,
    iBase_SILENT = iBase_ErrorActions_MIN,
    iBase_WARN_ONLY,
    iBase_THROW_ERROR,
    iBase_ErrorActions_MAX = iBase_THROW_ERROR
  };

    /*==========================================================
     * ERRORTYPE ENUMERATION
     *==========================================================
     */
  enum iBase_ErrorType {
    iBase_ErrorType_MIN = 0,
    iBase_SUCCESS = iBase_ErrorType_MIN,
    iBase_MESH_ALREADY_LOADED,
    iBase_NO_MESH_DATA,
    iBase_FILE_NOT_FOUND,
    iBase_FILE_WRITE_ERROR,
    iBase_NIL_ARRAY,
    iBase_BAD_ARRAY_SIZE,
    iBase_BAD_ARRAY_DIMENSION,
    iBase_INVALID_ENTITY_HANDLE,
    iBase_INVALID_ENTITY_COUNT,
    iBase_INVALID_ENTITY_TYPE,
    iBase_INVALID_ENTITY_TOPOLOGY,
    iBase_BAD_TYPE_AND_TOPO,
    iBase_ENTITY_CREATION_ERROR,
    iBase_INVALID_TAG_HANDLE,
    iBase_TAG_NOT_FOUND,
    iBase_TAG_ALREADY_EXISTS,
    iBase_TAG_IN_USE,
    iBase_INVALID_ENTITYSET_HANDLE,
    iBase_INVALID_ITERATOR_HANDLE,
    iBase_INVALID_ARGUMENT,
    iBase_MEMORY_ALLOCATION_FAILED,
    iBase_NOT_SUPPORTED,
    iBase_FAILURE,
    iBase_ErrorType_MAX = iBase_FAILURE
  };

    /*==========================================================
     * ERROR STRUCT
     *==========================================================
     */
  struct iBase_Error
  {
    int error_type;
    char description[120];
  };

    /*==========================================================
     * STORAGEORDER ENUMERATION
     *==========================================================
     */
  enum iBase_StorageOrder {
    iBase_StorageOrder_MIN = 0,
    iBase_BLOCKED = iBase_StorageOrder_MIN,
    iBase_INTERLEAVED,
    iBase_StorageOrder_MAX = iBase_INTERLEAVED
  };

    /*==========================================================
     * TAGVALUETYPE ENUMERATION
     *==========================================================
     */
  enum iBase_TagValueType {
    iBase_TagValueType_MIN = 0,
    iBase_BYTES = iBase_TagValueType_MIN,
    iBase_INTEGER,
    iBase_DOUBLE,
    iBase_ENTITY_HANDLE,
    iBase_ENTITY_SET_HANDLE,
    iBase_TagValueType_MAX = iBase_ENTITY_SET_HANDLE
  };

#ifdef __cplusplus
}
#endif

#endif /* #ifndef _ITAPS_iBase */
