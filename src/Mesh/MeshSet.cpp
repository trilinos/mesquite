// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD:  7-Nov-02 at 13:07:53 by Thomas Leurent
//
/*! \file MeshSet.cpp

\brief This files implements all the memory management issues related
to the copy of the original TSTT (or other maybe) mesh entity handles
into Mesquite.
That copy is of course encapsulated in the MeshSet class.
  
    \author Thomas Leurent
    \date 2002-05-16  
 */

#include "TSTT_Base.h"

#include "stdio.h"

#include "MeshSet.hpp"
#include "QualityImprover.hpp"


using namespace Mesquite;  


MeshSet::MeshSet() :
  spaceDim(0),
  mType(UNDEFINED_PATCH_TYPE),
  mParam1(0), mParam2(0),
  cullingMethodBits(0),
  fixedVertexTagName("fixed")
{
  currentMesh = meshSet.begin();
  currentVertex = verticesSet.begin();
}

// I have a feeling the destructor will do more later, so
// I didn't put it inline
MeshSet::~MeshSet()
{}


/*! \fn  MeshSet::add_mesh(TSTT::Mesh_Handle mh, MsqError &err)

    adds a TSTT mesh to the MeshSet. If used several times,
    it concatenates the vertices from several TSTT meshes into
    one Mesquite::MeshSet.
  */
#undef __FUNC__
#define __FUNC__ "MeshSet::add_mesh"
void MeshSet::add_mesh(TSTT::Mesh_Handle mh, MsqError &err)
{
  // sets MeshSet::SpaceDim
  int dim;
  TSTT::MeshError info=0;
  TSTT::Mesh_GetGeometricDimension(mh, &dim, &info);
  if (spaceDim == 0) 
    spaceDim = dim;
  else if (dim != spaceDim)
  {
    err.set_msg("Meshes of different dimensions added to the same MeshSet.");
    return;
  }
  
  // adds the TSTT::Mesh_Handle to the MeshSet.
  meshSet.push_front(mh);
}


/*! \fn MeshSet::reset(MsqError &err)

    Resets the MeshSet object.
    The current vertex is set back to the first vertex in the first mesh handle.
  */
#undef __FUNC__
#define __FUNC__ "MeshSet::reset"
void MeshSet::reset(MsqError &err)
{
  verticesSet.clear();
  currentVertex = verticesSet.begin();
  currentMesh = meshSet.begin();
}

bool MeshSet::set_patch_type(MeshSet::PatchType patch_type,
                             int patch_param1,
                             int patch_param2)
{
  mType = patch_type;
  mParam1 = patch_param1;
  mParam2 = patch_param2;
  
    // Setting the patch type restarts iterations through the mesh
  MsqError err;
  reset(err); MSQ_CHKERR(err);  
  
    // For now, we only support ELEMENTS_ON_VERTEX_PATCH
  if ( patch_type != ELEMENTS_ON_VERTEX_PATCH
       && patch_type != GLOBAL_PATCH )
  {
    return false;
  }
  return true;
}




/*! \fn MeshSet::get_next_vertices_set(MsqError &err)

    This function is usually called by MeshSet::get_next_patch().
    It retrieves the next vertices set from the TSTT interface.
    If the end has been reached, it resets the MeshSet and returns false .

    This function should be where memory management is optimized,
    by avoiding to fully duplicate the TSTT information.  

 */
#undef __FUNC__
#define __FUNC__ "MeshSet::get_next_vertices_set"
bool MeshSet::get_next_vertices_set(MsqError &err)
{
  // finds last mesh handle.
  std::list<TSTT::Mesh_Handle>::iterator last_mesh;
  last_mesh = meshSet.end();
  --last_mesh; // we want this to be the last Mesh handle, not one past the end.

  if (meshSet.empty()) {
    err.set_msg("MeshSet does not contain any mesh.");
    return false;
  }
 
  // if MeshSet is in initial state
  if ( verticesSet.empty() )
    currentMesh = meshSet.begin();
  // If we're at the end of the MeshSet::verticesSet but there is more
  else if ( currentVertex==verticesSet.end() && currentMesh!=last_mesh ) {
      currentMesh++;
  }
  else {
    // If we're at the end of the MeshSet
    if ( currentVertex==verticesSet.end() && currentMesh==last_mesh ) {
      // reset MeshSet
      currentMesh = meshSet.begin();
      verticesSet.clear();
      currentVertex = verticesSet.begin();
      return false; 
    }
  }
  
  // we need to create the next vertices set
  int num_vertices=0;
  TSTT::MeshError tstt_err=0;
  TSTT::Entity_Handle *vtx;
  TSTT::Mesh_GetEntities(*currentMesh, TSTT::VERTEX, &vtx, &num_vertices, &tstt_err);
  assert(!tstt_err);        
  /* Copy the vertices handles of one MeshHandle in the MeshSet::verticesSet.
     This can be modified to copy less entities
     and run with smaller merory footprint. */
  for (int i=0; i<num_vertices; ++i)
    {
      verticesSet.push_back(EntityEntry(*currentMesh, vtx[i]));
    }
  currentVertex = verticesSet.begin();
  // Can we free these here?  Are the handles useable in later
  // functions if we free them here?
  TSTT::Mesh_FreeEntityHandles(*currentMesh, vtx, &tstt_err);
  assert(!tstt_err);

  return true;
}


/*! \fn MeshSet::get_next_patch()
  \brief This function fills up a PatchData object with the mesh information
  necessary for optimization algorythms.
  The return value is true as long there exist a next patch, false otherwise.
  The list culling is performed in this function.

  This function is a friend of the PatchData class. Therefore the PatchData
  object passed as an argument is filled up directly.

  \param PatchData &pd: this is the PatchData object that will be filled up.
  \param int num_layers: number of layers of adjacencies included in the
  PatchData object. Set to -1 to retrieve the whole mesh.
*/
#undef __FUNC__
#define __FUNC__ "MeshSet::get_next_patch" 
bool MeshSet::get_next_patch(PatchData &pd,
                              MsqError &err )
{

  TSTT::MeshError tstt_err=0;

  // *************************************************
  // get the next vertices set from the TSTT interface.
  // This could also be replaced by a TSTT iterator.
  // *************************************************

  
  bool more_vtx;
  // if MeshSet object in initial state
  if (mType == ELEMENTS_ON_VERTEX_PATCH) {
    if ( verticesSet.empty() ) {
      more_vtx = get_next_vertices_set(err); MSQ_CHKERR(err);
      if ( more_vtx == false ) {
        return false; 
      }
    }
    else {
      currentVertex++;
      // If this is the end of the current vertices list
      if (currentVertex == verticesSet.end() ) {
        more_vtx = get_next_vertices_set(err);  MSQ_CHKERR(err);
        if ( more_vtx == false )
          return false; // no error needed, we're just at the end of the list.
      }
    }
  }
 
  MSQ_DEBUG_ACTION(3,{std::cout << "  o Patch Type: " << mType << std::endl; });
  
  // ***********************************************
  // Gathers the local patch information
  // ***********************************************

  if ( mType==ELEMENTS_ON_VERTEX_PATCH ) {

    // checks second argument.
    int num_layers = mParam1;
    if (num_layers != 1)
      {
        err.set_msg("no implementation for patch depth !=1 yet."); 
        return false;
      }
  
    const int MAX_NUM_VERTICES_PER_ELEM=12;
    TSTT::Entity_Handle* vertices = new TSTT::Entity_Handle[MAX_NUM_VERTICES_PER_ELEM];
    TSTT::Entity_Handle *patch_elements;
    TSTT::EntityTopology *element_topologies;
    int* csr_pointer; // = new int[500]; //dbg
    int* csr_data; // = new int[500]; //dbg
    //  int num_vertices=0;
    int num_elements;
    Vector3D vertex_coordV;
  
 //   ::set_new_handler(out_of_store);
  
    // This loop is the list culling process
    // also checks if this is the end of the vertices list
    // The culling assumes that you've got the same Mesh_Handle for all vertices
    // in verticesSet, so that the culling tag_handle stays the same.
    bool cull_vertex = true;
    // retrieves value of fixedVertexTagName .
    const char* bnd_tag_name;
    bnd_tag_name = fixedVertexTagName.c_str(); // converts from std::string to C-style string.
    void* bnd_tag_handle;
    if (cullingMethodBits & QualityImprover::NO_BOUNDARY_VTX) {
      TSTT::Mesh_tagGetHandle (currentVertex->mesh, bnd_tag_name, &bnd_tag_handle, &tstt_err);
    }
    while (cullingMethodBits!=0 && cull_vertex) {

      cull_vertex = false; // don't cull , by default

      // Culling of boudary vertices 
      if (cullingMethodBits & QualityImprover::NO_BOUNDARY_VTX) {
        int* on_boundary = NULL;
        int tag_size = sizeof(int);
        TSTT::Mesh_GetTag_Entity(currentVertex->mesh,
                                 (TSTT::cEntity_Handle) (currentVertex->entity),
                                 bnd_tag_handle, (void**)&on_boundary,
                                 &tag_size, &tstt_err);
      
        // Make sure the call succeeded.
        // NOTE: If the tag doesn't exist, we'll get an error...
        assert(!tstt_err);
        assert(on_boundary != NULL);
      
        if ( *on_boundary ==1)
          {
            MSQ_DEBUG_ACTION(2,{
              std::cout << "      o Culling vertex "
                        << distance(verticesSet.begin(), currentVertex)
                        << " according to NO_BOUNDARY_VTX." << std::endl; });
            cull_vertex = true;
            currentVertex++;
            // If this is the end of the MeshSet vertices list
            if (currentVertex == verticesSet.end() ) {
              more_vtx = get_next_vertices_set(err);  MSQ_CHKERR(err);
              if ( more_vtx == false )
                return false; // no error needed, we're just at the end of the list.
            }
          }
      
        // Delete memory allocated by Mesh_GetTag_Entity.
        delete on_boundary;
      }
      // other type of culling ...
      //else if () {
      //}

    } // end of culling loop
  


      // retrieves array of pointers to regions adjacent to the free vertex ...
    num_elements = 0; // this will make TSTT allocate the patch_regions array
    TSTT::Entity_GetAdjacencies(currentVertex->mesh,
                                (TSTT::cEntity_Handle *)&(currentVertex->entity),
                                1, TSTT::REGION, &patch_elements,
                                &csr_pointer, &csr_data, 
                                &num_elements, &tstt_err );
    // ... if there are no regions adjacent to the vtx,
    // retrieves adjacent faces.

    if ( num_elements == 0 )
      {
        TSTT::Entity_GetAdjacencies(currentVertex->mesh,
                                    (TSTT::cEntity_Handle *) &(currentVertex->entity),
                                    1, TSTT::FACE, &patch_elements,
                                    &csr_pointer, &csr_data,
                                    &num_elements, &tstt_err);
      }
    assert(!tstt_err);
    //std::cout << "num_elements: " << num_elements << std::endl;
    // If there are neither regions or faces adjacent to the vertex, EXIT_FAILURE.
    if ( num_elements == 0 )
      {
        err.set_msg("no regions or faces are adjacent to free vertex.");
        return false;
      }
  
    // Now, get the types of the adjacent elements.
    element_topologies = new enum TSTT::EntityTopology[num_elements];
    TSTT::Entity_GetTopology( currentVertex->mesh,
                              (TSTT::cEntity_Handle *) patch_elements, num_elements,
                              element_topologies, &tstt_err);
  
    // finds an upper bound to the number of vertices in the patch.
    // The worst case (2 triangles in the patch with center vertex
    // on the boundary) must work.
    int max_num_vertices = 1; // center vtx
    for (int i=0; i<num_elements; ++i)
      {
        int n = MsqMeshEntity::vertex_count(element_topologies[i], err); MSQ_CHKERR(err);
        max_num_vertices += n-1; // no need to count center vertex several times.
      }

    // gets rid of previous Patch information (but keeps memory allocated).
    pd.clear();
    // allocates memory for coordinates array
    // if current memory is inapropriate
    pd.reserve_vertex_capacity(max_num_vertices, err); MSQ_CHKERR(err);
  
    // allocates memory for connectivity array
    // if current memory is inapropriate
    pd.reserve_element_capacity(num_elements, err); MSQ_CHKERR(err);
  
    int num_coords=3*MAX_NUM_VERTICES_PER_ELEM;
    double* vertex_coords = new double[num_coords];
    // enters the coordinates of the central vertex in PatchData::coordsArray
    // at position 0.
    TSTT::Entity_GetVertexCoords(currentVertex->mesh,
                                 (TSTT::cEntity_Handle *) &(currentVertex->entity),
                                 1, TSTT::INTERLEAVED,
                                 &vertex_coords, &num_coords, &tstt_err);
    pd.add_vertex(currentVertex->mesh, currentVertex->entity,
                  &vertex_coords[0], false, err); MSQ_CHKERR(err);
    pd.set_num_free_vertices(1);
  
    // For each element within the patch,
    // fill up the adequate info in PatchData.
    for (int e=0; e<num_elements; ++e)
      {
        // gets the vertices of an element ...
        int num_element_vtx=MAX_NUM_VERTICES_PER_ELEM;
        TSTT::Entity_GetAdjacencies( currentVertex->mesh,
                                     (TSTT::cEntity_Handle *) &(patch_elements[e]),
                                     1, TSTT::VERTEX, &vertices,
                                     &csr_pointer, &csr_data, 
                                     &num_element_vtx, &tstt_err );
        TSTT::Entity_GetVertexCoords(currentVertex->mesh,
                                     (TSTT::cEntity_Handle *) vertices,
                                     num_element_vtx, TSTT::INTERLEAVED,
                                     &vertex_coords, &num_coords,
                                     &tstt_err);
        assert(!tstt_err);
        // ... enters the coordinates of those vertices in PatchData ...
        int vtx_ind[MAX_NUM_VERTICES_PER_ELEM];
        for (int n=0; n<num_element_vtx; ++n)
          {
            vtx_ind[n] = pd.add_vertex(currentVertex->mesh, vertices[n],
                                       &vertex_coords[3*n], true, err); MSQ_CHKERR(err);
          }
        // ... and adds the element to PatchData
        pd.add_element(currentVertex->mesh, patch_elements[e],
                       vtx_ind, tstt_to_mesquite(element_topologies[e]), err); MSQ_CHKERR(err);
      }
    delete[] element_topologies;
    delete[] vertices;
    delete[] vertex_coords;
  
    //cout << "numVertices for local patch: " << pd.num_vertices() << endl;
  
    // TODO : provide patch for several layers of adjacencies.
  
    // free entities allocated in TSTT::Entity_GetAdjacencies()

    TSTT::Mesh_FreeEntityHandles(currentVertex->mesh, patch_elements, &tstt_err);
  
    return true;
  }    

    // ***********************************************
    // Gathers the global patch information
    // ***********************************************
  
  else if (mType==GLOBAL_PATCH) {
  
    const int MAX_NUM_VERTICES_PER_ELEM=12;
    TSTT::Entity_Handle* elem_vtx = new TSTT::Entity_Handle[MAX_NUM_VERTICES_PER_ELEM];
    TSTT::Entity_Handle *all_elements;
    TSTT::EntityTopology *element_topologies;
    int* csr_pointer; 
    int* csr_data;
    int num_elements;
    int num_vertices;
//    Vector3D vertex_coordV;
  
 //   ::set_new_handler(out_of_store);

    // loops over all TSTT mesh handles
    for ( currentMesh = meshSet.begin();
          currentMesh != meshSet.end();
          ++currentMesh ) {
      
      // retrieves all the elements of the mesh
      num_elements = 0; // this will make TSTT allocate the patch_regions array
      TSTT::Mesh_GetEntities(*currentMesh,
                             TSTT::REGION, &all_elements,
                             &num_elements, &tstt_err );
      // ... if there are no regions adjacent to the vtx,
      // retrieves adjacent faces.
      if ( num_elements == 0 )
        {
          TSTT::Mesh_GetEntities(*currentMesh,
                                 TSTT::FACE, &all_elements,
                                 &num_elements, &tstt_err );
        }
      assert(!tstt_err);

      // If there are neither regions or faces adjacent to the vertex, EXIT_FAILURE.
      if ( num_elements == 0 ) {
        err.set_msg("no regions or faces are adjacent to free vertex.");
        return false;
      }
  
      // Now, get the types of the elements.
      element_topologies = new enum TSTT::EntityTopology[num_elements];
      TSTT::Entity_GetTopology( *currentMesh,
                                (TSTT::cEntity_Handle *) all_elements, num_elements,
                                element_topologies, &tstt_err);
  
      // finds the number of vertices in the patch
      TSTT::Mesh_Services_GetInt( *currentMesh, "number of vertices",
                            &num_vertices, &tstt_err);

      // gets rid of previous Patch information (but keeps memory allocated).
      pd.clear();
      // allocates memory for coordinates array
      // if current memory is inapropriate
      pd.reserve_vertex_capacity(num_vertices, err); MSQ_CHKERR(err);
  
      // allocates memory for connectivity array
      // if current memory is inapropriate
      pd.reserve_element_capacity(num_elements, err); MSQ_CHKERR(err);
  
      int elem_num_coords=3*MAX_NUM_VERTICES_PER_ELEM;
      double* elem_vtx_coords = new double[elem_num_coords];


      enum store_in_patch {
        free_vtx,
        fixed_vtx,
        over
      };
      
      // retrieves value of fixedVertexTagName.
      // The value can be different for each TSTT mesh handle.
      const char* bnd_tag_name;
      bnd_tag_name = fixedVertexTagName.c_str(); // converts from std::string to C-style string.
      void* bnd_tag_handle;
      TSTT::Mesh_tagGetHandle(*currentMesh, bnd_tag_name, &bnd_tag_handle, &tstt_err);
      
      // For each element within the patch,
      // fill up the adequate info in PatchData.
      store_in_patch store = free_vtx;
      while (store != over) {
        // before starting 2nd pass
        if (store == fixed_vtx) {
          int num_free_vtx = pd.num_vertices(); // gets nb of vtx added in 1st pass
          pd.set_num_free_vertices(num_free_vtx);
        }
        
        for (int e=0; e<num_elements; ++e) {
          
          // gets the vertices of an element ...
          int num_element_vtx=MAX_NUM_VERTICES_PER_ELEM;
          TSTT::Entity_GetAdjacencies( *currentMesh,
                                       (TSTT::cEntity_Handle *) &(all_elements[e]),
                                       1, TSTT::VERTEX, &elem_vtx,
                                       &csr_pointer, &csr_data, 
                                       &num_element_vtx, &tstt_err );
          TSTT::Entity_GetVertexCoords(*currentMesh,
                                       (TSTT::cEntity_Handle *) elem_vtx,
                                       num_element_vtx, TSTT::INTERLEAVED,
                                       &elem_vtx_coords, &elem_num_coords,
                                       &tstt_err);
          assert(!tstt_err);
          // ... enters the coordinates of those vertices in PatchData ...
          int vtx_ind[MAX_NUM_VERTICES_PER_ELEM];
          for (int n=0; n<num_element_vtx; ++n)
            {
              // retrieves the fixed tad value
              int* on_boundary = NULL;
              int tag_size = sizeof(int);
              TSTT::Mesh_GetTag_Entity(*currentMesh,
                                       (TSTT::cEntity_Handle) elem_vtx[n],
                                       bnd_tag_handle, (void**)&on_boundary,
                                       &tag_size, &tstt_err);
              assert(!tstt_err);
              assert(on_boundary != NULL);
              if (   store==free_vtx && *on_boundary==0    // first pass we add free vertices
                  || store==fixed_vtx ) { // second pass we add all vertices
                                        // (this retrieves index of previously added vertices). 
                vtx_ind[n] = pd.add_vertex(*currentMesh, elem_vtx[n],
                                           &elem_vtx_coords[3*n], true, err); MSQ_CHKERR(err);
              }
            }
          // ... and adds the element to PatchData if in second pass.
          if (store==fixed_vtx) {
            pd.add_element(*currentMesh, all_elements[e],
                           vtx_ind, tstt_to_mesquite(element_topologies[e]), err);
            MSQ_CHKERR(err);
          }
        }
        
        if (store == fixed_vtx)  store=over; // for next round
        else if (store == free_vtx)   store=fixed_vtx; // for next round
      }

      delete[] element_topologies;
      delete[] elem_vtx;
      delete[] elem_vtx_coords;
  
      //cout << "numVertices for local patch: " << pd.num_vertices() << endl;
  
      // TODO : provide patch for several layers of adjacencies.
  
      // free entities allocated in TSTT::Entity_GetAdjacencies()
      TSTT::Mesh_FreeEntityHandles(*currentMesh, all_elements, &tstt_err);
    }
  
    return true;
  }  

  
  err.set_msg("no implementation for specified patch type.");
  return false;

}

#include "MeanRatioQualityMetric.hpp"
#include "LInfTemplate.hpp"
#include "ConjugateGradient.hpp"

bool MeshSet::get_next_element_group(PatchData &pd, MsqError &err)
{
    //VERY temp solution just to test
    //ShapeQualityMetric* sm = MeanRatioQualityMetric::create_new();
    //LInfTemplate li=LInfTemplate(sm);
    //ConjugateGradient cg=ConjugateGradient(&li);
    //cg.set_patch_depth(1);
   return get_next_patch(pd, err);
     //delete sm;
     //delete &li;
     //delete &cg;
}

bool MeshSet::get_next_node_group(PatchData &pd, MsqError &err)
{

  err.set_msg("no implementation yet."); MSQ_CHKERR(err);
  
  return false; // ... no implementation 
  
}


// ************* AOMD tmp TEST *************

// struct PrintVertex {
//   void operator()(cEntity_Ptr const m)
//   {
//     int entityID=1,info=1;
//     Entity_GetID(m,&entityID,&info);
//     std::cout << "Vertex " << entityID << "\t" << info << std::endl;
//   }
// };
 
// void Mesquite::test_aomd(void)
// {
//   Mesh* theMesh = newMesh();
 
//   int info;
//   // Load the mesh
//   MeshLoad(theMesh,"cube.sms",&info);
//   Iter vIter = beginIter(theMesh,0);
//   Iter vEnd  = endIter(theMesh,0);
 
//   std::for_each(vIter,vEnd,PrintVertex());
//   MeshSave(theMesh,"cube.msh",&info);
//   return;
// }
