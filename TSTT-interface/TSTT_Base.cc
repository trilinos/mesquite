/* ****************************************************************************
 *
 *  /TSTT/TSTT_Base.cc
 *  Created by Eunyoung Seol on Tue Jul 10 2002, 14:04:31 EDT
 *  Modified by Yunhua Luo on Aug 7 - Aug 12, 2002
 *  
 *  Changes in this version mainly include: 
 *    (1) places where the specifications are not fulfilled;
 *    (2) the results from the teleconference on Aug. 8, 2002.
 *
 *  SCOREC
 *
 *
 *  File Contents: implementation file of TSTT_Base.h
 *
 *  Version: %I%
 *  Date modified: %G%
 *
 *************************************************************************** */

#include "TSTT_Base.h"

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <stdlib.h>
#include <assert.h>
 
#include "mMesh.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mPoint.h"
#include "mAttachableDataContainer.h"
#include "mAOMD.h"
#include "AOMD.h"
#include "AOMDInternals.h"
 
using AOMD::AOMD_Util;
using AOMD::mMesh;
using AOMD::mEntity;
using AOMD::mPoint;
using AOMD::mVertex;
using std::string;
using std::cout;
using std::endl;
 
namespace TSTT{
const char* _0D = "number of 0D entities"; 
const char* _1D = "number of 1D entities"; 
const char* _2D = "number of 2D entities"; 
const char* _3D = "number of 3D entities"; 
const char* vertices = "number of vertices"; 
const char* edges = "number of edges"; 
const char* faces = "number of faces";    
const char* regions = "number of regions"; 

//===================================================================
void  Mesh_Create(Mesh_Handle* mesh, MeshError *errorReturn)
//===================================================================
{
  *errorReturn = (void*)1;
  *mesh = new mMesh;
  *errorReturn = NULL;
}

//===================================================================
void  Mesh_Load(cMesh_Handle mesh, const char *name, MeshError *errorReturn)
//===================================================================
{
  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;
  AOMD_Util::Instance()->import(name, m);
  Int dim = (m->size(3))? 3 : ((m->size(2)) ? 2 : 1);
  m->modifyState(dim, dim-1, true, 0);
  if (dim == 3)
    m->modifyState(2,1,true, 0);
  *errorReturn = NULL;
}

//===================================================================
void  Mesh_Services_GetInt(cMesh_Handle mesh, const char* key, 
//===================================================================
                    Int *value, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;

  if (!strcmp(key, _0D) || !strcmp(key, vertices))
    *value = m->size(0);
  else if (!strcmp(key, _1D) || !strcmp(key, edges))
    *value = m->size(1);
  else if (!strcmp(key, _2D) || !strcmp(key, faces))
    *value = m->size(2);
  else if (!strcmp(key, _3D) || !strcmp(key, regions))
    *value = m->size(3);
  else
    *value = -1;
  *errorReturn = NULL;
}

//===================================================================
void  Mesh_Services_GetDouble(cMesh_Handle mesh, const char* key, 
//===================================================================
                    Double *value, MeshError *errorReturn)
{ 
  *errorReturn = (void*) 1;
  //mMesh* m = (mMesh*) mesh;

  //if (!strcmp(key, _0D) || !strcmp(key, vertices))
  //  *value = m->size(0);
  //else if (!strcmp(key, _1D) || !strcmp(key, edges))
  //  *value = m->size(1);
  //else if (!strcmp(key, _2D) || !strcmp(key, faces))
  //  *value = m->size(2);
  //else if (!strcmp(key, _3D) || !strcmp(key, regions))
  //  *value = m->size(3);
  //else
  //  *value = -1.0;
   cout << "No string key values for this function are supported at this time." << endl;
  *errorReturn = NULL;
}

//===================================================================
void  Mesh_Services_GetFloat(cMesh_Handle mesh, const char* key,  
//===================================================================
                    Float *value, MeshError *errorReturn)
{  
  *errorReturn = (void*)1;
  //mMesh* m = (mMesh*) mesh;   
  // No string key values for this function are supported at this time. (yluo) 
  //if (!strcmp(key, _0D) || !strcmp(key, vertices))
  //  *value = m->size(0);
  //else if (!strcmp(key, _1D) || !strcmp(key, edges))
  //  *value = m->size(1);
  //else if (!strcmp(key, _2D) || !strcmp(key, faces))
  //  *value = m->size(2);
  //else if (!strcmp(key, _3D) || !strcmp(key, regions))
  //  *value = m->size(3);
  //else
  //  *value = -1.0;
  cout << "No string key values for this function are supported at this time." << endl;
  *errorReturn = NULL;
}

//===================================================================
void Mesh_Services_GetString(cMesh_Handle mesh, const char* key, 
//===================================================================
                    char** value, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  //mMesh* m = (mMesh*) mesh;
  cout << "No string key values for this function are supported at this time." << endl;
  *errorReturn = NULL;
}

//===================================================================
void Mesh_GetGeometricDimension(cMesh_Handle mesh, Int *dim, 
//===================================================================
                                   MeshError *errorReturn)
{ 
  // "Note that this (dim) does not necessarily equal to the highest topological
  //  entity as surface meshes consisting of two dimensional topological 
  //  elements can live in three-dimension space." (from TSTT_Base.h). Yes
  //  this is correct, but with the above assumption, the geometric dimension
  //  can only be determined by the model which uses it, a mesh database can 
  //  not tell how it would be used by a model. The following is a 
  //  temporary implementation, but also might be the only way to implement
  //  this function. (yluo)
  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;
  *dim = (m->size(3))? 3 : ((m->size(2)) ? 2 : 1 );
  *errorReturn = NULL; 
}

//===================================================================
void  Mesh_GetEntities(cMesh_Handle mesh, const EntityType entity_type,
//=================================================================== 
          Entity_Handle **entity_handles, Int *num_entities, 
          MeshError *errorReturn)
{ 
  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;

  Int numEnt = m->size(entity_type);
  if (numEnt == 0) return;                                 // (yluo)
  if (*num_entities > 0 && *num_entities < numEnt) return; // (yluo)
  if (*num_entities == 0)  // The mesh allocates the array
    *entity_handles = new Entity_Handle[numEnt];
  //(yluo)else if (numEnt > 0 && *num_entities < numEnt)
  //(yluo)return;
  // else enough memory has been allocated to *entity_handles already 
  
  // the mesh fills the array
  mMesh::iter it = m->begin(entity_type);
  mMesh::iter itend = m->end(entity_type);
  Int counter = 0;
  for (; it!=itend; ++it)
    (*entity_handles)[counter++] = (void*)*it;

  *num_entities = numEnt;

  *errorReturn = NULL; 
}   
   
//=================================================================== 
void  Mesh_GetEntityAdjacency(cMesh_Handle mesh, 
//===================================================================
                const EntityType entity_requestor,
                const EntityType adjacency_entity_requested,
                Entity_Handle **entity_handles, 
                Int **csr_pointer, Int **csr_data,
                Int *num_adj_entities,
                MeshError *errorReturn)
{
  *errorReturn = (void*)1; 
  mMesh* m = (mMesh*) mesh;
  
  if ((int)entity_requestor < (int)adjacency_entity_requested)
     m->modifyState((int)entity_requestor, (int)adjacency_entity_requested,true,0);
  mMesh::iter it = m->begin(entity_requestor);
  mMesh::iter itend = m->end(entity_requestor);
  if (*it == 0) 
  {  cout<<"FATAL: the entity_requestor doesn't exist.\n";
    return;
  }

  // This should be changed later when higher-order adjacency
  // is available. (yluo) 
  if (entity_requestor == adjacency_entity_requested)
  { 
    cout<<"FATAL: entity_requestor should be different from adjacency_entity_requested\n";
    return;
  }
 
  Int counter=0, nb=0;
  Int numAdjEntities=0;
  for (;it!=itend;++it)
    numAdjEntities += (*it)->size(adjacency_entity_requested);

  //cout<<"numAdjEntities = "<<numAdjEntities<<endl;
  // Do memory allocation
  if (numAdjEntities <= 0) return;  // (yluo)
  if (*num_adj_entities != 0 && *num_adj_entities < numAdjEntities) return; // (yluo)
  if (*num_adj_entities==0)
    *entity_handles = new Entity_Handle[numAdjEntities];
  //(yluo) else if (*num_adj_entities < numAdjEntities)
  //(yluo)   return;
  
  *csr_pointer = new Int[m->size(entity_requestor)];
  *csr_data = new Int[m->size(entity_requestor)]; 
  *num_adj_entities = numAdjEntities;

  it = m->begin(entity_requestor);
  for (;it!=itend; ++it)
  {
    Int numAdj = (*it)->size(adjacency_entity_requested);
    (*csr_pointer)[nb] = numAdj;
    (*csr_data)[nb] = counter;
    for (Int i=0; i<numAdj; ++i,++counter)
      (*entity_handles)[counter] = (*it)->get(adjacency_entity_requested,i);
    ++nb;
  }
  assert(counter == numAdjEntities);
  assert(nb == m->size(entity_requestor)); 
  *errorReturn = NULL;
}

// WHY DO WE NEED MESH HANDLE FOR THIS FUNCTION?
// WE CAN GET PREFERRED STORAGE ORDER WITHOUT MESH HANDLE!!!
//===================================================================
void Mesh_GetNativeStorageOrder(cMesh_Handle mesh,
//=================================================================== 
                   StorageOrder *storage_order, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  //mMesh* m = (mMesh*) mesh;
  *storage_order = TSTT::UNDETERMINED;
  *errorReturn = NULL;
}

//===================================================================
void Mesh_GetCoordinates(cMesh_Handle mesh, Int* num_entities, 
//===================================================================
                         Double** coordinates, StorageOrder *storage_order,
                         char *options, MeshError *errorReturn)
{ 
  // Based on the agreement achieved at the teleconference on Aug. 8 (Thu),
  // 2002, "range(n1:n2)", "interleaved" and "blocked" were removed from 
  // "option". So the implementation of this function was re-written.
  // The options "local" and "includeHigherOrderNodes" are not implemented
  // at this time for reasons from both definition and implementation 
  // technique (yluo)

  int i=0, j=0;
  mMesh* m = (mMesh*) mesh;
  double ent_coords[3];
  int nV = m->size(VERTEX);
  
  if (nV == 0) return;
  if (*num_entities != 0 && *num_entities < nV) return;
  if (*num_entities == 0)
  { 
     *coordinates = new Double[nV*3];
     *num_entities = nV;
  }

  mMesh::iter it = m->begin(VERTEX);
  mMesh::iter itend = m->end(VERTEX);
  
  for( ; it != itend; ++it )
  {
    mEntity *e = *it;
    mPoint *p = ((mVertex*)e)->ppoint();
    ent_coords[0] = (*p)(0);
    ent_coords[1] = (*p)(1);
    ent_coords[2] = (*p)(2);

    if (*storage_order == INTERLEAVED)
    {
       for (i=0; i<3; i++)
         (*coordinates)[j++] = ent_coords[i];
    }
    else
    {
      (*coordinates)[j] = ent_coords[0];    	
      (*coordinates)[nV + j] = ent_coords[1];    	
      (*coordinates)[2*nV + j] = ent_coords[2];    	
       j++;
    }
  }
  *errorReturn = (void*)NULL;
}

//===================================================================
void Mesh_GetEntityIndicies(cMesh_Handle mesh, 
//===================================================================
                   const EntityType entityType, 
                   Int** count, Int* num_count,
                   Int** index, Int* num_index, 
                   char *options, MeshError *errorReturn )
{
  // The two arguments "num_count" and "num_index" were added here
  // purely for C language. For SIDL version, these two variables
  // are not necessary. The implementation of this function has been
  // changed so that if "entityType"="VERTEX", the function returns
  // the global IDs of all vertices, according the result of 
  // discussion on e-note (page 5, answers from Lori). This function
  // was originally designed to be used with other functions, e.g.
  // Mesh_GetCoordinates(...) to get the IDs of vertices. (yluo)
  //
  // Another point is that, although "expanded" option are also given 
  // in the specification, we here chose "csr" as it might be the only
  // feasible way for mixed meshes where the numbers of element nodes
  // are different. (yluo)

  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;

  mMesh::iter it = m->begin(entityType);
  mMesh::iter itend = m->end(entityType);
  if (*it==NULL)
    return;
  Int num_vertices=0, num_entities=0;
  for (; it!=itend; ++it)
  { 
    mEntity* e = *it;
    if (entityType == VERTEX)
      num_vertices++;
    else
      num_vertices += e->size(VERTEX);
    ++num_entities;
  }

  if (num_entities == 0 || num_vertices ==0) return;
  if (*num_count != 0 && *num_count < num_entities) return;
  if (*num_index != 0 && *num_index < num_vertices) return;
 
  // Do Memory Allocation IF NEEDED
  if (*num_count == 0)
  {
    *count = new Int[num_entities];
    *num_count = num_entities;
  }
  if (*num_index==0)
  {  
    *index = new Int[num_vertices];
    *num_index = num_vertices;
  }

  // Parse input options
  bool csr = true;      // WE ASSUME CSR IS DEFAULT FORMAT
  EntityTopology topo = TSTT::UNDEFINED;
/*
  bool includeHigherOrderNodes = false;
  string sOptions(options);

  // The following two options are implemented at this time
  // for both definition and implementation reasions. (yluo) 
  if (sOptions.find("expanded") != string::npos) 
    csr = false;
  if (sOptions.find("includeHigherOrderNodes") != string::npos)
    includeHigherOrderNodes = true;

  // The following implementation seems better than the options 
  // described in the specification. With the following implementation, 
  // users can inquiry elements with any supported topology type, rather
  // than only "HEXAHEDRON". But this option value does not work for the
  // reasons described later. (yluo)
    if (sOptions.find("POLYGON")!=string::npos)
      topo = TSTT::POLYGON;   
    else if (sOptions.find("TRIANGLE")!=string::npos) 
      topo = TSTT::TRIANGLE;
    else if (sOptions.find("QUADRILATERAL")!=string::npos)
      topo = TSTT::QUADRILATERAL;
    else if (sOptions.find("POLYHEDRON")!=string::npos)
      cout<<"ALERT: AOMD doesn't support POLYHEDRON\n";  
    else if (sOptions.find("TETRAHEDRON")!=string::npos)
      topo = TSTT::TETRAHEDRON;
    else if (sOptions.find("HEXAHEDRON")!=string::npos)
      topo = TSTT::HEXAHEDRON;
    else if (sOptions.find("PRISM")!=string::npos)
      topo = TSTT::PRISM;
    else if (sOptions.find("PYRAMID")!=string::npos)
      topo = TSTT::PYRAMID;
    else if (sOptions.find("SEPTAHEDRON")!=string::npos)
      cout<<"ALERT: AOMD doesn't support SEPTAHEDRON\n";  
*/
  if (csr)
  {   
    it = m->begin(entityType);
    Int pos_count=0, pos_index=0;

    if (topo == TSTT::UNDEFINED)
    {
      for (; it!=itend; ++it)
      {
        mEntity* e = *it;
        (*count)[pos_count++] = pos_index;
        Int numVertex;
        if (entityType == VERTEX) numVertex = 1;
        else numVertex = e->size(TSTT::VERTEX);
        for (Int j=0; j<numVertex; j++)
        {
          mEntity* v;
          if (entityType == VERTEX) v = *it;
          else v = e->get(TSTT::VERTEX, j);
          (*index)[pos_index++] = v->getId();
        }
        (*count)[pos_count] = pos_index;
      }  
      assert(pos_count==*num_count && pos_index == *num_index);
    } //TSTT::UNDEFINED
    else
    {
      /*for (; it!=itend; ++it)
      {
        mEntity* e = *it;
        // There are quite several problems to be fixed here:
        // (1) the definitions of "EntityType" and "TopologyType"
        //     in TSTT_Base.h are inconsistent with those in AOMD.
        // (2) function mEntity::getType() does not work. It needs
        //     quite several changes in AOMD to fix this problem.
        // for the above reasons, especially the second one, the 
        // option part is not implemented this time. (yluo)
        if ((Int)((*it)->getType()) == (Int)topo)
        {
          cout<<"topo = "<<topo<<endl;
          (*count)[pos_count++] = pos_index;
          Int numVertex;
          if (entityType == VERTEX)
	    numVertex = 1;
          else
            numVertex = e->size(TSTT::VERTEX);
          for (Int j=0; j<numVertex; ++j)
          {
            mEntity* v;
            if (entityType == VERTEX) v = *it;
            else v = e->get(TSTT::VERTEX, j);
            (*index)[pos_index++] = v->getId();
          }
          (*count)[pos_count] = pos_index;
        }
      }
      *num_count = pos_count;
      *num_index = pos_index;
      */
      *num_count = 0;
      *num_index = 0;
      cout << "option = 'TRI, TET or HEX, ...' is not implemented yet." << endl;
      } // !TSTT::UNDEFINED
  }  // CSR
  else // expanded
  {

  } 
  *errorReturn = NULL;
}
  
//===================================================================
void Mesh_GetEntityAdjacencyIndicies(cMesh_Handle mesh, 
//===================================================================
               const EntityType entityType, 
               const EntityType entityAdjacencyType,
               Int **count, Int* num_count,
               Int **index, Int* num_index,
               char *options, MeshError *errorReturn)
{
  // The two arguments "num_count" and "num_index" were added purely for
  // C language, as there is no way to know the length of an array in C
  // except explicitly given. (yluo) 
  //
  // option = "csr" is assumed, the other options are not implemented 
  // this time. (yluo) 

  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;
  mMesh::iter it = m->begin(entityType);
  mMesh::iter itend = m->end(entityType);
  if (entityType < entityAdjacencyType)
      m->modifyState(entityType,entityAdjacencyType,true,0);
  Int numAdjEnt=0;
  Int numEnt = m->size(entityType);
  for (; it!=itend; ++it)
    numAdjEnt += (*it)->size(entityAdjacencyType);
 
  if (numAdjEnt == 0 || numEnt == 0) return;
  if (*num_count != 0 && *num_count < numEnt) return;
  if (*num_index != 0 && *num_index < numAdjEnt) return;

  // Do Memory Allocation IF NEEDED
  if (*num_count == 0)  { *count = new Int[numEnt];    *num_count = numEnt; }
  if (*num_index==0)    { *index = new Int[numAdjEnt]; *num_index = numAdjEnt; }

  // Parse input options
  bool csr = true;   	// WE ASSUME CSR IS DEFAULT FORMAT
  bool includeHigherOrderNodes = false;

  /* 
  string sOptions(options); 
  if (sOptions.find("expanded") != string::npos)
    csr = false;
  if (sOptions.find("includeHigherOrderNodes") != string::npos)
    includeHigherOrderNodes = true;
  */

  it = m->begin(entityType);

  if (csr) 
  {  
    Int pos_index = 0;
    Int pos_count = 0;
    for (; it!=itend; ++it)
    { 
      Int numAdj = (*it)->size(entityAdjacencyType);
      (*count)[pos_count++] = pos_index;
      for (Int i=0; i<numAdj; ++i)
       (*index)[pos_index++] = (*it)->get(entityAdjacencyType, i)->getId();
    }
    (*count)[pos_count] = pos_index;
    assert(pos_count==*num_count && pos_index==*num_index);
  }  // CSR
  else // expanded
  {
  // IT WILL BE DONE LATER
  }
  *errorReturn = NULL;
}

// 2 FUNCTIONS REGARDING WORKSET ARE NOT IMPLEMENTED AT THIS TIME, 
// MORE DISCUSSION WITH YUNHUA IS NEEDED
//===================================================================
//void Mesh_InitializeWorksetIterator(cMesh_Handle mesh,
//===================================================================
//              const EntityType entity_type,
//              const Int requested_workset_size,
//              Entity_Workset_Iterator *workset_iterator,
//              MeshError *errorReturn)
//{
//  *errorReturn = (void*)1;
//  mMesh* m = (mMesh*) mesh;
//  *errorReturn = NULL;
//}

//===================================================================
//void Mesh_GetNextWorkset(cMesh_Handle mesh, Entity_Workset_Iterator workset_iterator,
//===================================================================
//           Entity_Handle **entity_handles, Int *num_entities,
//           MeshError *errorReturn)
//{
//  *errorReturn = (void*)1;
//  *errorReturn = NULL;
//}

// "WorkSet" is not defined in AOMD. As a temporary solution, here we provide
//  the following two functions.
//====================================================================
void Mesh_InitializeEntityIterator(cMesh_Handle mesh, 
//====================================================================
             EntityType entityType, 
	     Entity_Handle *firstEntity,
	     Entity_Handle *lastEntity, 
	     MeshError *errorReturn)
{
	*errorReturn = (void*)1;
	mMesh* m = (mMesh*) mesh;
	mMesh::iter it = m->begin(entityType);
        mMesh::iter itend = m->end(entityType);
	*firstEntity = (void*)(*it);
	*lastEntity = (void*)(*itend);
	*errorReturn = NULL;
}

//======================================================================
void Mesh_GetNextEntity(cMesh_Handle mesh, EntityType entityType,
//======================================================================
           Entity_Handle *entity, MeshError *errorReturn)
{
	*errorReturn = (void*)1;
	mMesh* m = (mMesh*) mesh;
	mMesh::iter it = m->begin(entityType);
        mMesh::iter itend = m->end(entityType);
	for (; it != itend; ++it)
	{
		if ( *it == *entity )
		{
			it++;
			*entity = (void*)(*it);
			break;
		}
	}
	*errorReturn = NULL;
}
//===================================================================
void  Mesh_FreeEntityHandles(cMesh_Handle mesh, Entity_Handle *entity_handles,
//===================================================================
                MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  delete[] entity_handles;
  *errorReturn = NULL;
}

//===================================================================
void  Mesh_FreeEntityAdjacency(cMesh_Handle mesh, Entity_Handle *entity_handles, 
//===================================================================
                Int *csr_pointer, Int *csr_data, 
                MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  delete [] entity_handles;
  delete [] csr_pointer;
  delete [] csr_data;
  *errorReturn = NULL;
}

//===================================================================
void  Mesh_Destroy(Mesh_Handle mesh, MeshError *errorReturn)
//===================================================================
{
  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;
  delete mesh;
  mesh = (mMesh*)0;
  *errorReturn = NULL; 
}

//===================================================================
void  Mesh_GetEntitiesFromID(cMesh_Handle mesh, Int* ID, 
//===================================================================
	          EntityType entity_type, const Int num_identifiers,	  
              Entity_Handle **entity_handle, 
              MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  mMesh* m = (mMesh*) mesh;
  for (Int i=0; i<num_identifiers; ++i)
  {
    mMesh::iter it = m->begin(entity_type); 
    mMesh::iter itend = m->end(entity_type);
    for (; it!=itend; ++it)
    {
      if ((Int)ID[i] == (*it)->getId())
      {
        (*entity_handle)[i] =(void*)*it;
        break;
      }
    }
  }
  *errorReturn = NULL;
}

//===================================================================
void  Entity_GetGlobalID(cMesh_Handle mesh, cEntity_Handle *entity_handles, 
//===================================================================
                Int num_entities, Int **ID, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  mEntity** e = (mEntity**) entity_handles;
  for (Int i=0; i<num_entities; ++i, ++e)
     (*ID)[i] = (*e)->getId();
  *errorReturn = NULL;
}

/*
This function has gone based on the discussion on Aug. 8, 2002. (yluo)
//===================================================================
void  Entity_GetGeometricDimension(cEntity_Handle *entity_handles, 
//===================================================================
                const Int num_entities,
                Int *dim, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  mEntity** e = (mEntity**)entity_handles; 
  for (Int i=0; i<num_entities; ++i, ++e)
    dim[i] = (*e)->getClassification()->getLevel();
  *errorReturn = NULL;

}
*/

//===================================================================
void  Entity_GetTopologicalDimension(cMesh_Handle mesh, cEntity_Handle *entity_handles, 
//===================================================================
                const Int num_entities,
                Int *dim, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  mEntity** e = (mEntity**)entity_handles;
  // RETURNS ONE OF 0,1,2,3
  for (Int i=0; i<num_entities; ++i, ++e)
    dim[i] = (*e)->getLevel();
  *errorReturn = NULL;
}
    
// ACTUALLY, WE BELIEVE THAT EntityType and EntityTopology SHOULD BE CHANGED
// IN AOMD, TOPOLOGY AND DIMENSION ARE IDENTICAL
//===================================================================
void  Entity_GetTopology(cMesh_Handle mesh, cEntity_Handle *entity_handles, 
//===================================================================
                const Int num_entities,
                EntityTopology *topology,
                MeshError *errorReturn)
{
  // this function dose not work normally for the following reason:
  //   (1) the definition of "EntityTyep" and "EntityTopology" in 
  //       TSTT is inconsistent with that in AOMD. (this is relatively
  //       easier to solve)
  //   (2) mEntity::getType() is not well designed in AOMD and does 
  //       not work correctly. (yluo)
  *errorReturn = (void*)1;
  mEntity** e = (mEntity**) entity_handles;
  // RETURNS ONE OF VERTEX,EDGE,FACE,TRI,QUAD,TET,HEX,PRISM,PYRAMID
  for (Int i=0; i<num_entities; ++i, ++e)
    {
      AOMD::mEntity::mType topo = (*e)->getType();
      switch(topo) {
      case AOMD::mEntity::VERTEX:
        topology[i] = TSTT::POINT;
        break;
      case AOMD::mEntity::EDGE:
        topology[i] = TSTT::LINE;
        break;
      case AOMD::mEntity::TRI: 
        topology[i] = TSTT::TRIANGLE;
        break;
      case AOMD::mEntity::QUAD:
        topology[i] = TSTT::QUADRILATERAL;
        break;
      case AOMD::mEntity::HEX:
        topology[i] = TSTT::HEXAHEDRON;
        break;
      case AOMD::mEntity::PRISM:
        topology[i] = TSTT::PRISM;
        break;
      case AOMD::mEntity::PYRAMID:
        topology[i] = TSTT::PYRAMID;
        break;
      case AOMD::mEntity::TET:
        topology[i] = TSTT::TETRAHEDRON;
        break;
      default:
        return;
      }
    }
  *errorReturn = NULL;
}
 
//===================================================================
void Entity_GetType(cMesh_Handle mesh, cEntity_Handle *entity_handles, 
//===================================================================
               const Int num_entities, 
               EntityType* type, 
               MeshError* errorReturn)
{
  // it is not clear that what is the essential difference between 
  // this function and function Entity_GetToplogicalDimension(...)
  // (yluo)
  *errorReturn = (void*)1;
  mEntity** e = (mEntity**) entity_handles;
  // RETURNS ONE OF VERTEX,EDGE,FACE,REGION
  for (Int i=0; i<num_entities; ++i, ++e)
    type[i] = (EntityType)(*e)->getLevel();
  *errorReturn = NULL;
}

//===================================================================
void  Entity_GetVertexCoords(cMesh_Handle mesh, cEntity_Handle *entity_handles,
//=================================================================== 
	           Int num_entities, StorageOrder storage_order,
               Double **coords, Int *num_coords, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  if ( num_entities == 0) return;
  if (*num_coords != 0 && *num_coords < num_entities*3)  return;
  if (*num_coords == 0) // allocate memory as needed
  {   *coords = new Double[num_entities*3]; 
      *num_coords = num_entities*3;
  }
  mEntity** e = (mEntity**) entity_handles;
 
  Int count=0;
  if (storage_order == TSTT::BLOCKED)    
    for (Int i=0; i<num_entities;++i,++e) 
    {
      mPoint* p = ((mVertex*)*e)->ppoint();
      (*coords)[count] = (*p)(0);
      (*coords)[num_entities+count] = (*p)(1);
      (*coords)[2 * num_entities+count] =(*p)(2);
      ++count;
    } // if (blocked)
  else // (storage_order==TSTT::INTERLEAVED)
    for (Int i=0; i<num_entities; ++i, ++e)
    {
      mPoint* p = ((mVertex*) *e)->ppoint();
      for (Int k=0; k<3; ++k)
        (*coords)[count++] = (*p)(k);
    }
  *errorReturn = NULL;
}


//===================================================================
void  Entity_SetVertexCoords(cMesh_Handle mesh, cEntity_Handle *entity_handles,
//===================================================================
                             Int num_entities, StorageOrder storage_order, Int num_dimension,
                             Double *coords, MeshError *errorReturn)
{
  *errorReturn = (void*)1;
  if ( num_entities == 0) return;
  if ( num_dimension!=2 && num_dimension!=3 ) return;

  mEntity** e = (mEntity**) entity_handles;

  Int count=0;
  Double x, y, z;
  if (storage_order == TSTT::BLOCKED)
    for (Int i=0; i<num_entities;++i,++e)
    {
      x = coords[count];
      y = coords[num_entities+count];
      if (num_dimension==3)
        z = coords[2 * num_entities+count];
      else // ==2
        z = 0;
      mPoint p(x,y,z);
      ((mVertex*)*e)->move(p); // moves vertex to new position p
      ++count;
    } // if (blocked)
  else // (storage_order==TSTT::INTERLEAVED)
    for (Int i=0; i<num_entities; ++i, ++e)
    {
      x = coords[count++];
      y = coords[count++];
      if (num_dimension==3)
        z = coords[count++];
      else // ==2
        z = 0;
      mPoint p(x,y,z);
      ((mVertex*) *e)->move(p);
    }
  *errorReturn = NULL;
}


//===================================================================
void  Entity_GetAdjacencies(cMesh_Handle mesh, cEntity_Handle *entity_handles,
//===================================================================
               Int num_entities,
               const EntityType entity_type, 
               Entity_Handle **adj_entities_handle, 
 	       Int **csr_pointer, Int **csr_data,
               Int *num_adjacent,
               MeshError *errorReturn)
{ 
  *errorReturn = (void*)1;
  if (num_entities == 0) return;
  mEntity** e = (mEntity**) entity_handles;
  mMesh *m = (mMesh*)mesh;
  int requestorType, requestedType;
  requestorType = (*e)->getLevel();
  requestedType = (int)entity_type;
  if (requestorType < requestedType)
     m->modifyState(requestorType, requestedType,true,0);
  Int numAdjEntities=0;
  for (Int i=0;i<num_entities;++i,++e)
    numAdjEntities += (*e)->size(entity_type);
  if (numAdjEntities == 0) return;
  if (*num_adjacent != 0 && *num_adjacent < numAdjEntities)
    return;
  // Do memory allocation
  if (*num_adjacent==0)
  {
    *adj_entities_handle = new Entity_Handle[numAdjEntities];
  }
  *num_adjacent = numAdjEntities;
  
  *csr_pointer = new Int[num_entities];
  *csr_data = new Int[num_entities]; 
  
  Int counter=0, nb=0;
  e = (mEntity**) entity_handles;
  for (Int i=0; i<num_entities; ++i, ++e)
  {
    Int numAdj = (*e)->size(entity_type);
    (*csr_pointer)[nb] = numAdj;
    (*csr_data)[nb] = counter;
    for (Int j=0; j<numAdj; ++j,++counter)
      (*adj_entities_handle)[counter] = (*e)->get(entity_type,j);
    ++nb;
  }
  *errorReturn = NULL;
}


//======================================================
  void Mesh_tagGetHandle (const char *tag_name,
//======================================================
                     void** tag_handle,
                     MeshError *errorReturn)
{
  pMeshDataId* tagId = new pMeshDataId;
  *tagId = MD_lookupMeshDataId(tag_name);
  *tag_handle = (void*) tagId;
}


//======================================================
 void Mesh_GetTag_Entity(cMesh_Handle mesh,
//======================================================                         
                         cEntity_Handle this_entity,
                         void* tag_handle,
                         void **tag_value, Int *tag_size,
                         MeshError *errorReturn)
{
  *errorReturn = (void*)0;
  mEntity* e = (mEntity*) this_entity;

   int* value = new int;
   pMeshDataId* tagId = (pMeshDataId*) tag_handle;
   EN_getDataInt(e, *tagId, value);
  
   *tag_size = sizeof(int);
   *tag_value = (void*) value;
}


}  // end of namespace TSTT

