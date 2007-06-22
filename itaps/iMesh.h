#ifndef IMESH_CBIND_H__
#define IMESH_CBIND_H__

#ifndef ITAPS
#define ITAPS
#endif

#include "iBase.h"
#include "iMesh_protos.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef void* iMesh_Instance;
  typedef void* iMesh_EntityIterator;
  typedef void* iMesh_EntityArrIterator;

  enum iMesh_EntityTopology {
    iMesh_POINT = 0,              /**< a general zero-dimensional entity  */
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
    iMesh_ALL_TOPOLOGIES      /**< allows the user to request information
			       *     about all the topology types */
  };

  void iMesh_getErrorType(iMesh_Instance instance,
                          /*out*/ int *error_type, 
                          int *err);

  void iMesh_getDescription(iMesh_Instance instance,
                            /*inout*/ char *descr, 
                            int *err, 
                            /*in*/ int descr_len);

  void iMesh_newMesh(const char *options,
                     /*out*/ iMesh_Instance *instance, 
                     /*out*/ int *err, 
                     /*in*/ int options_len);

  void iMesh_dtor(iMesh_Instance instance, 
                  /*out*/ int *err);

  void iMesh_load(iMesh_Instance instance,
                  /*in*/ const iBase_EntitySetHandle entity_set_handle,
                  /*in*/ const char *name, 
                  /*in*/ const char *options,
                  /*out*/ int *err, 
                  /*in*/ int name_len, 
                  /*in*/ int options_len);

  void iMesh_save(iMesh_Instance instance,
                  /*in*/ const iBase_EntitySetHandle entity_set_handle,
                  /*in*/ const char *name, 
                  /*in*/ const char *options,
                  /*out*/ int *err, 
                  /*in*/ const int name_len, 
                  /*in*/ int options_len);

  void iMesh_getRootSet(iMesh_Instance instance,
                        /*out*/ iBase_EntitySetHandle *root_set, 
                        /*out*/ int *err);

  void iMesh_getGeometricDimension(iMesh_Instance instance,
                                   /*out*/ int *geom_dim, 
                                   /*out*/ int *err);

  void iMesh_getDfltStorage(iMesh_Instance instance,
                            /*out*/ int *order, 
                            /*out*/ int *err);

  void iMesh_getAdjTable (iMesh_Instance instance,
                          /*out*/ int** adjacency_table,
                          /*inout*/ int* adjacency_table_allocated,
                          /*out*/ int* adjacency_table_size, 
                          /*out*/ int *err);

  void iMesh_getNumOfType(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_type,
                          /*out*/ int *num_type, 
                          /*out*/ int *err);

  void iMesh_getNumOfTopo(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_topology,
                          /*out*/ int *num_topo, 
                          /*out*/ int *err);

  void iMesh_areEHValid(iMesh_Instance instance, 
                        /*in*/ int doReset,
                        /*out*/ int *areHandlesInvariant, 
                        /*out*/ int *err);

  void iMesh_getAllVtxCoords (iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set_handle,
                              /*inout*/ double** coordinates,
                              /*inout*/ int* coordinates_allocated,
                              /*out*/ int* coordinates_size,
                              /*inout*/ int** in_entity_set,
                              /*inout*/ int* in_entity_set_allocated,
                              /*out*/ int* in_entity_set_size,
                              /*inout*/ int* storage_order, /*out*/ int *err);

  void iMesh_getVtxCoordIndex (iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle entity_set_handle,
                               /*in*/ const int requested_entity_type,
                               /*in*/ const int requested_entity_topology,
                               /*in*/ const int entity_adjacency_type,
                               /*inout*/ int** offset,
                               /*inout*/ int* offset_allocated,
                               /*out*/ int* offset_size,
                               /*inout*/ int** index,
                               /*inout*/ int* index_allocated,
                               /*out*/ int* index_size,
                               /*inout*/  int** entity_topologies,
                               /*inout*/ int* entity_topologies_allocated,
                               /*out*/ int* entity_topologies_size, 
                               /*out*/ int *err);

  void iMesh_getEntities(iMesh_Instance instance,
                         /*in*/ const iBase_EntitySetHandle entity_set_handle,
                         /*in*/ const int entity_type,
                         /*in*/ const int entity_topology,
                         /*out*/ iBase_EntityHandle** entity_handles,
                         /*out*/ int* entity_handles_allocated,
                         /*out*/ int* entity_handles_size,
                         /*out*/ int *err);

  void iMesh_getVtxArrCoords(iMesh_Instance instance,
                             /*in*/ const iBase_EntityHandle* vertex_handles,
                             /*in*/ const int vertex_handles_size,
                             /*inout*/ int* storage_order,
                             /*inout*/ double** coords,
                             /*inout*/ int* coords_allocated,
                             /*out*/ int* coords_size, 
                             /*out*/ int *err);


  void iMesh_getAdjEntities(iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle entity_set_handle,
                            /*in*/ const int entity_type_requestor,
                            /*in*/ const int entity_topology_requestor,
                            /*in*/ const int entity_type_requested,
                            /*inout*/ iBase_EntityHandle** adj_entity_handles,
                            /*inout*/ int* adj_entity_handles_allocated,
                            /*out*/ int* adj_entity_handles_size,
                            /*inout*/ int** offset,
                            /*inout*/ int* offset_allocated,
                            /*out*/ int* offset_size,
                            /*inout*/ int** in_entity_set,
                            /*inout*/ int* in_entity_set_allocated,
                            /*out*/ int* in_entity_set_size, 
                            /*out*/ int *err);

  void iMesh_initEntArrIter(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const int requested_entity_type,
                            /*in*/ const int requested_entity_topology,
                            /*in*/ const int requested_array_size,
                            /*out*/ iMesh_EntityArrIterator* entArr_iterator,
                            /*out*/ int *err);

  void iMesh_getNextEntArrIter(iMesh_Instance instance,
                               /*in*/ iMesh_EntityArrIterator entArr_iterator,
                               /*inout*/ iBase_EntityHandle** entity_handles,
                               /*inout*/ int* entity_handles_allocated,
                               /*out*/ int* entity_handles_size,
                               /*out*/ int *has_data, 
                               /*out*/ int *err);


  void iMesh_resetEntArrIter(iMesh_Instance instance,
                             /*in*/ iMesh_EntityArrIterator entArr_iterator, 
                             /*out*/ int *err);


  void iMesh_endEntArrIter(iMesh_Instance instance,
                           /*in*/ iMesh_EntityArrIterator entArr_iterator, 
                           /*out*/ int *err);


  void iMesh_getEntArrTopo(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** topology,
                           /*inout*/ int* topology_allocated,
                           /*out*/ int* topology_size, 
                           /*out*/ int *err);


  void iMesh_getEntArrType(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** type,
                           /*inout*/ int* type_allocated,
                           /*out*/ int* type_size, 
                           /*out*/ int *err);


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

  void iMesh_createEntSet(iMesh_Instance instance,
                          /*in*/ const int isList,
                          /*out*/ iBase_EntitySetHandle* entity_set_created,
                          /*out*/ int *err);


  void iMesh_destroyEntSet(iMesh_Instance instance,
                           /*in*/ iBase_EntitySetHandle entity_set,
                           /*out*/ int *err);

  void iMesh_isList(iMesh_Instance instance,
                    /*in*/ const iBase_EntitySetHandle entity_set,
                    /*out*/ int *is_list,
                    /*out*/ int *err);

  void iMesh_getNumEntSets(iMesh_Instance instance,
                           /*in*/ const iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const int num_hops,
                           /*out*/ int *num_sets,
                           /*out*/ int *err);


  void iMesh_getEntSets(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set_handle,
                        /*in*/ const int num_hops,
                        /*out*/ iBase_EntitySetHandle** contained_set_handles,
                        /*out*/ int* contained_set_handles_allocated,
                        /*out*/ int* contained_set_handles_size,
                        /*out*/ int *err);

  void iMesh_addEntToSet(iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle entity_handle,
                         /*inout*/ iBase_EntitySetHandle* entity_set,
                         /*out*/ int *err);

  void iMesh_rmvEntFromSet(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle entity_handle,
                           /*inout*/ iBase_EntitySetHandle* entity_set,
                           /*out*/ int *err);


  void iMesh_addEntArrToSet(iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ const int entity_handles_size,
                            /*inout*/ iBase_EntitySetHandle* entity_set,
                            /*out*/ int *err);


  void iMesh_rmvEntArrFromSet(iMesh_Instance instance,
                              /*in*/ const iBase_EntityHandle* entity_handles,
                              /*in*/ const int entity_handles_size,
                              /*inout*/ iBase_EntitySetHandle* entity_set,
                              /*out*/ int *err);


  void iMesh_addEntSet(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_set_to_add,
                       /*inout*/ iBase_EntitySetHandle* entity_set_handle,
                       /*out*/ int *err);


  void iMesh_rmvEntSet(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle entity_set_to_remove,
                       /*inout*/ iBase_EntitySetHandle* entity_set_handle,
                       /*out*/ int *err);

  void iMesh_isEntContained(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle containing_entity_set,
                            /*in*/ const iBase_EntitySetHandle contained_entity,
                            /*out*/ int *is_contained,
                            /*out*/ int *err);

  void iMesh_isEntSetContained(iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle containing_entity_set,
                               /*in*/ const iBase_EntitySetHandle contained_entity_set,
                               /*out*/ int *is_contained,
                               /*out*/ int *err);

  void iMesh_addPrntChld(iMesh_Instance instance,
                         /*inout*/ iBase_EntitySetHandle* parent_entity_set,
                         /*inout*/ iBase_EntitySetHandle* child_entity_set,
                         /*out*/ int *err);

  void iMesh_rmvPrntChld(iMesh_Instance instance,
                         /*inout*/ iBase_EntitySetHandle* parent_entity_set,
                         /*inout*/ iBase_EntitySetHandle* child_entity_set,
                         /*out*/ int *err);

  void iMesh_isChildOf(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle parent_entity_set,
                       /*in*/ const iBase_EntitySetHandle child_entity_set,
                       /*out*/ int *is_child,
                       /*out*/ int *err);

  void iMesh_getNumChld(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        /*out*/ int *num_child,
                        /*out*/ int *err);

  void iMesh_getNumPrnt(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        /*out*/ int *num_parent,
                        /*out*/ int *err);

  void iMesh_getChldn(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*out*/ iBase_EntitySetHandle** entity_set_handles,
                      /*out*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size,
                      /*out*/ int *err);

  void iMesh_getPrnts(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*out*/ iBase_EntitySetHandle** entity_set_handles,
                      /*out*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size,
                      /*out*/ int *err);

  void iMesh_setVtxArrCoords(iMesh_Instance instance,
                             /*in*/ iBase_EntityHandle* vertex_handles,
                             /*in*/ const int vertex_handles_size,
                             /*in*/ const int storage_order,
                             /*in*/ const double* new_coords,
                             /*in*/ const int new_coords_size,
                             /*out*/ int *err);


  void iMesh_createVtxArr(iMesh_Instance instance,
                          /*in*/ const int num_verts,
                          /*in*/ const int storage_order,
                          /*in*/ const double* new_coords,
                          /*in*/ const int new_coords_size,
                          /*inout*/ iBase_EntityHandle** new_vertex_handles,
                          /*inout*/ int* new_vertex_handles_allocated,
                          /*inout*/ int* new_vertex_handles_size,
                          /*out*/ int *err);


  void iMesh_createEntArr(iMesh_Instance instance,
                          /*in*/ const int new_entity_topology,
                          /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                          /*in*/ const int lower_order_entity_handles_size,
                          /*out*/ iBase_EntityHandle** new_entity_handles,
                          /*out*/ int* new_entity_handles_allocated,
                          /*out*/ int* new_entity_handles_size,
                          /*inout*/ int** status,
                          /*inout*/ int* status_allocated,
                          /*out*/ int* status_size,
                          /*out*/ int *err);


  void iMesh_deleteEntArr(iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*out*/ int *err);


  void iMesh_createTag(iMesh_Instance instance,
                       /*in*/ const char* tag_name,
                       /*in*/ const int tag_size,
                       /*in*/ const int tag_type,
                       /*out*/ iBase_TagHandle* tag_handle, 
                       /*out*/ int *err,
                       /*in*/ const int tag_name_len);


  void iMesh_destroyTag(iMesh_Instance instance,
                        /*in*/ iBase_TagHandle tag_handle,
                        /*in*/ const int forced,
                        /*out*/ int *err);

  void iMesh_getTagName(iMesh_Instance instance,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*inout*/ char *name,
                        /*out*/ int *err,
                        /*in*/ int name_len);

  void iMesh_getTagSizeValues(iMesh_Instance instance,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ int *tag_size,
                              /*out*/ int *err);

  void iMesh_getTagSizeBytes(iMesh_Instance instance,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*out*/ int *tag_size,
                             /*out*/ int *err);

  void iMesh_getTagHandle(iMesh_Instance instance,
                          /*in*/ const char* tag_name,
                          /*out*/ iBase_TagHandle *tag_handle, 
                          /*out*/ int *err,
                          int tag_name_len);

  void iMesh_getTagType(iMesh_Instance instance,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ int *tag_type,
                        /*out*/ int *err);


  void iMesh_setEntSetData(iMesh_Instance instance,
                           /*in*/ iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const char* tag_value,
                           /*in*/ const int tag_value_size,
                           /*out*/ int *err);


  void iMesh_setEntSetIntData(iMesh_Instance instance,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*in*/ const int tag_value,
                              /*out*/ int *err);


  void iMesh_setEntSetDblData(iMesh_Instance instance,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*in*/ const double tag_value,
                              /*out*/ int *err);


  void iMesh_setEntSetEHData(iMesh_Instance instance,
                             /*in*/ iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*in*/ const iBase_EntityHandle tag_value,
                             /*out*/ int *err);


  void iMesh_getEntSetData(iMesh_Instance instance,
                           /*in*/ const iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ char** tag_value,
                           /*inout*/ int* tag_value_allocated,
                           /*inout*/ int* tag_value_size,
                           /*out*/ int *err);

  void iMesh_getEntSetIntData(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ int *out_data,
                              /*out*/ int *err);

  void iMesh_getEntSetDblData(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ double *out_data,
                              /*out*/ int *err);

  void iMesh_getEntSetEHData(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*out*/ iBase_EntityHandle *out_data,
                             /*out*/ int *err);

  void iMesh_getAllEntSetTags(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set_handle,
                              /*out*/ iBase_TagHandle** tag_handles,
                              /*out*/ int* tag_handles_allocated,
                              /*out*/ int* tag_handles_size,
                              /*out*/ int *err);

  void iMesh_rmvEntSetTag(iMesh_Instance instance,
                          /*in*/ iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*out*/ int *err);

  void iMesh_setVtxCoords(iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle vertex_handle,
                          /*in*/ const double x, /*in*/ const double y,
                          /*in*/ const double z,
                          /*out*/ int *err);

  void iMesh_createVtx(iMesh_Instance instance,
                       /*in*/ const double x, /*in*/ const double y,
                       /*in*/ const double z,
                       /*out*/ iBase_EntityHandle* new_vertex_handle,
                       /*out*/ int *err);

  void iMesh_createEnt(iMesh_Instance instance,
                       /*in*/ const int new_entity_topology,
                       /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                       /*in*/ const int lower_order_entity_handles_size,
                       /*out*/ iBase_EntityHandle* new_entity_handle,
                       /*out*/ int* status,
                       /*out*/ int *err);

  void iMesh_deleteEnt(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       /*out*/ int *err);

  void iMesh_getArrData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle* entity_handles,
                        /*in*/ const int entity_handles_size,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*inout*/ char** tag_values,
                        /*inout*/int* tag_values_allocated,
                        /*out*/ int* tag_values_size,
                        /*out*/ int *err);

  void iMesh_getIntArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ int** tag_values,
                           /*inout*/ int* tag_values_allocated,
                           /*out*/ int* tag_values_size,
                           /*out*/ int *err);

  void iMesh_getDblArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ double** tag_values,
                           /*inout*/ int* tag_values_allocated,
                           /*out*/ int* tag_values_size,
                           /*out*/ int *err);

  void iMesh_getEHArrData(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*inout*/ iBase_EntityHandle** tag_value,
                          /*inout*/ int* tag_value_allocated,
                          /*out*/ int* tag_value_size,
                          /*out*/ int *err);

  void iMesh_setArrData(iMesh_Instance instance,
                        /*in*/ iBase_EntityHandle* entity_handles,
                        /*in*/ const int entity_handles_size,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const char* tag_values,
                        /*in*/ const int tag_values_size,
                        /*out*/ int *err);

  void iMesh_setIntArrData(iMesh_Instance instance,
                           /*in*/ iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const int* tag_values,
                           /*in*/ const int tag_values_size,
                           /*out*/ int *err);

  void iMesh_setDblArrData(iMesh_Instance instance,
                           /*in*/ iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const double* tag_values,
                           /*in*/ const int tag_values_size,
                           /*out*/ int *err);

  void iMesh_setEHArrData(iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*in*/ const iBase_EntityHandle* tag_values,
                          /*in*/ const int tag_values_size,
                          /*out*/ int *err);

  void iMesh_rmvArrTag(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle* entity_handles,
                       /*in*/ const int entity_handles_size,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*out*/ int *err);

  void iMesh_getData(iMesh_Instance instance,
                     /*in*/ const iBase_EntityHandle entity_handle,
                     /*in*/ const iBase_TagHandle tag_handle,
                     /*inout*/ char** tag_value,
                     /*inout*/ int *tag_value_allocated,
                     /*out*/ int *tag_value_size,
                     /*out*/ int *err);

  void iMesh_getIntData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ int *out_data,
                        /*out*/ int *err);

  void iMesh_getDblData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ double *out_data,
                        /*out*/ int *err);

  void iMesh_getEHData(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*out*/ iBase_EntityHandle *out_data,
                       /*out*/ int *err);

  void iMesh_setData(iMesh_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ const iBase_TagHandle tag_handle,
                     /*in*/ const char* tag_value,
                     /*in*/ const int tag_value_size,
                     /*out*/ int *err);

  void iMesh_setIntData(iMesh_Instance instance,
                        /*in*/ iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const int tag_value,
                        /*out*/ int *err);

  void iMesh_setDblData(iMesh_Instance instance,

                        /*in*/ iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const double tag_value,
                        /*out*/ int *err);

  void iMesh_setEHData(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*in*/ const iBase_EntityHandle tag_value,
                       /*out*/ int *err);

  void iMesh_getAllTags(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*inout*/ iBase_TagHandle** tag_handles,
                        /*inout*/ int* tag_handles_allocated,
                        /*out*/ int* tag_handles_size,
                        /*out*/ int *err);

  void iMesh_rmvTag(iMesh_Instance instance,
                    /*in*/ iBase_EntityHandle entity_handle,
                    /*in*/ const iBase_TagHandle tag_handle,
                    /*out*/ int *err);

  void iMesh_initEntIter(iMesh_Instance instance,
                         /*in*/ const iBase_EntitySetHandle entity_set_handle,
                         /*in*/ const int requested_entity_type,
                         /*in*/ const int requested_entity_topology,
                         /*out*/ iMesh_EntityIterator* entity_iterator,
                         /*out*/ int *err);

  void iMesh_getNextEntIter(iMesh_Instance instance,
                            /*in*/ iMesh_EntityIterator entity_iterator,
                            /*out*/ iBase_EntityHandle* entity_handle,
                            /*out*/ int *has_data,
                            /*out*/ int *err);

  void iMesh_resetEntIter(iMesh_Instance instance,
                          /*in*/ iMesh_EntityIterator entity_iterator,
                          /*out*/ int *err);

  void iMesh_endEntIter(iMesh_Instance instance,
                        /*in*/ iMesh_EntityIterator entity_iterator,
                        /*out*/ int *err);

  void iMesh_getEntTopo(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*out*/ int *out_topo,
                        /*out*/ int *err);

  void iMesh_getEntType(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*out*/ int *out_type,
                        /*out*/ int *err);

  void iMesh_getVtxCoord(iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle vertex_handle,
                         /*out*/ double *x, /*out*/ double *y, /*out*/ double *z,
                         /*out*/ int *err);

  void iMesh_getEntAdj(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const int entity_type_requested,
                       /*inout*/ iBase_EntityHandle** adj_entity_handles,
                       /*inout*/ int* adj_entity_handles_allocated,
                       /*out*/ int* adj_entity_handles_size,
                       /*out*/ int *err);

  void iMesh_subtract(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle entity_set_1,
                      /*in*/ const iBase_EntitySetHandle entity_set_2,
                      /*out*/ iBase_EntitySetHandle* result_entity_set,
                      /*out*/ int *err);

  void iMesh_intersect(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle entity_set_1,
                       /*in*/ const iBase_EntitySetHandle entity_set_2,
                       /*out*/ iBase_EntitySetHandle* result_entity_set,
                       /*out*/ int *err);

  void iMesh_unite(iMesh_Instance instance,
                   /*in*/ const iBase_EntitySetHandle entity_set_1,
                   /*in*/ const iBase_EntitySetHandle entity_set_2,
                   /*out*/ iBase_EntitySetHandle* result_entity_set,
                   /*out*/ int *err);

#ifdef __cplusplus
}
#endif

#endif // ifndef IMESH_CBIND_H__
