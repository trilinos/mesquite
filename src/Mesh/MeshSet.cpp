// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD:  3-Feb-04 at 13:00:56 by Thomas Leurent
//
/*! \file MeshSet.cpp

\brief This files implements all the memory management issues related
to the copy of the original TSTT (or other maybe) mesh entity handles
into Mesquite.
That copy is of course encapsulated in the MeshSet class.
  
    \author Thomas Leurent
    \date 2002-05-16  
 */

#include "MeshSet.hpp"
#include "QualityImprover.hpp"

using namespace Mesquite;  

MSQ_USE(cout);
MSQ_USE(endl);

MeshSet::MeshSet() :
  vertexIterator(NULL),
  spaceDim(0),
  csrOffsets(0),
  csrData(0),
  vertArray(NULL),
  elemArray(NULL),
  elemTopologies(0),
  vertexOnBoundary(0),
  csrOffsetsSize(0),
  csrDataSize(0),
  vertArraySize(0),
  elemArraySize(0),
  elemTopologiesSize(0),
  mDomain(NULL)
{
  cullFlag=MsqVertex::MSQ_SOFT_FIXED;
}

MeshSet::~MeshSet()
{
    // Delete the vertex iterator
  delete vertexIterator;
    // Release all of our meshes
  list<Mesquite::Mesh*>::iterator it = meshSet.begin();
  while(!(it == meshSet.end()))
    (*it++)->release();
  
    // Delete cache arrays
  delete [] csrOffsets;
  delete [] csrData;
  delete [] vertArray;
  delete [] elemArray;
  delete [] elemTopologies;
  delete [] vertexOnBoundary;
}


/*! \fn  MeshSet::add_mesh(Mesquite::Mesh* mh, MsqError &err)

    adds a Mesquite::Mesh to the MeshSet. If used several times,
    it concatenates the vertices from several meshes into
    one Mesquite::MeshSet.
  */
#undef __FUNC__
#define __FUNC__ "MeshSet::add_mesh"
void MeshSet::add_mesh(Mesquite::Mesh* mesh, MsqError &err)
{
    // sets MeshSet::SpaceDim
  int dim = mesh->get_geometric_dimension(err); MSQ_CHKERR(err);
  if (spaceDim == 0) 
    spaceDim = dim;
  else if (dim != spaceDim)
  {
    err.set_msg("Meshes of different dimensions added to the same MeshSet.");
    return;
  }
  
    // adds the Mesh* to the MeshSet.
  meshSet.push_front(mesh);
}


/*! \fn MeshSet::reset(MsqError &err)

    Resets the MeshSet object.
    The current vertex is set back to the first vertex in the first mesh handle.
  */
#undef __FUNC__
#define __FUNC__ "MeshSet::reset"
void MeshSet::reset(MsqError& err)
{
    // If we have at least one mesh...
  if (meshSet.size())
  {
      // If we aren't already on the first mesh...
    if (!vertexIterator || !(currentMesh == meshSet.begin()))
    {
      currentMesh = meshSet.begin();
      delete vertexIterator;
      vertexIterator = (*currentMesh)->vertex_iterator(err); MSQ_CHKERR(err);
    }
    else // We ARE on the first mesh...
    {
        // ...so we can re-use the iterator.
      vertexIterator->restart();
    }
  }
  else
  {
      // This is probably redundant...
    vertexIterator = NULL;
  }
}


/*! \fn MeshSet::get_next_patch(PatchData &pd, MsqError &err )
  \brief This function fills up a PatchData object with the mesh information
  necessary for optimization algorithms.
  The return value is true as long there exists a next patch, false otherwise.
  The list culling is performed in this function.
  
  This function is a friend of the PatchData class. Therefore the PatchData
  object passed as an argument is filled up directly.
  
  \param PatchData &pd: this is the PatchData object that will be filled up.
*/
#undef __FUNC__
#define __FUNC__ "MeshSet::get_next_patch" 
bool MeshSet::get_next_patch(PatchData &pd,
                             PatchDataParameters &pd_params,
                             MsqError &err )
{
  FUNCTION_TIMER_START(__FUNC__);

    // get rid of previous Patch information (but keep memory allocated).
  pd.clear();

    // Mark this MeshSet as the originator
  pd.meshSet = this;
  if (mDomain != NULL)  pd.domainSet = true;
  else                  pd.domainSet = false;
    // Get the patch parameters.
  PatchData::PatchType patch_type = pd_params.get_patch_type();
  long unsigned int culling_method_bits = pd_params.get_culling_method_bits();
  
  MSQ_DEBUG_ACTION(3,{cout << "  o Patch Type: "
                                << patch_type << endl; });

  switch (patch_type)
  {
    case PatchData::ELEMENTS_ON_VERTEX_PATCH:
    {
        //variable to store the center vertex's fixed flag
      MsqVertex::FlagMask center_fixed_byte;
        // Make sure we're only getting a patch depth of 1
      int num_layers = pd_params.get_nb_layers(err); MSQ_CHKERR(err);
      if (num_layers != 1)
      {
        err.set_msg("no implementation for patch depth !=1."); 
        FUNCTION_TIMER_END();  
        return false;
      }

        // Set the patch type
      pd.mType = PatchData::ELEMENTS_ON_VERTEX_PATCH;
      
        // If this is our first time through the mesh,
        // initialize everything.
      if (!vertexIterator)
        reset(err);
      
        // currentVertex is pointing at next potential center vertex.
        // Move forward in the list of vertices if necessary.
      bool next_vertex_identified = false;
      while (!next_vertex_identified)
      {
          // Move to next mesh if necessary
        if (vertexIterator->is_at_end())
        {
          delete vertexIterator;
          ++currentMesh;
          if (currentMesh == meshSet.end())
          {
            vertexIterator = NULL;
            FUNCTION_TIMER_END();  
            return false;
          }
          vertexIterator = (*currentMesh)->vertex_iterator(err); MSQ_CHKERR(err);
        }
          //if this is a 'boundary' fixed flag, skip it now
        else if ((culling_method_bits & PatchData::NO_BOUNDARY_VTX))
        {
          bool on_bnd;
          Mesquite::Mesh::VertexHandle vtx = **vertexIterator;
          (*currentMesh)->vertices_are_on_boundary(&vtx, &on_bnd, 1, err);
          if (on_bnd)
            ++(*vertexIterator);
        }
          //otherwise we check to see if this vertex has been culled
        else{
            //get the fixed_bit_flag for the center vertex
          (*currentMesh)->vertex_get_byte(**vertexIterator,&center_fixed_byte, err);
            //remove the hard fixed flag if it has been set
          center_fixed_byte &= ~(MsqVertex::MSQ_HARD_FIXED);
            //if it is culled, skip it
          if(center_fixed_byte & cullFlag)
          {
            ++(*vertexIterator);
          }
          else
          {
              // We found the right one
            next_vertex_identified = true;
          }//end else (vertex was not fixed [boundary] or culled)
        }//end else (iterator was not at the end and vertex was not boundary)  
      }//end while (!next_vertex_identified)
      Mesh::VertexHandle vertex = **vertexIterator;
      ++(*vertexIterator);
      
        // Get the number of elements in this vertex
      size_t num_elems =
        (*currentMesh)->vertex_get_attached_element_count(vertex, err);
      MSQ_CHKERR(err);
      
        // Get the elements attached to this vertex
      if (elemArraySize < num_elems)
      {
        delete [] elemArray;
        elemArray = new Mesh::ElementHandle[num_elems];
        elemArraySize = num_elems;
      }
      
      (*currentMesh)->vertex_get_attached_elements(vertex,
                                                   elemArray,
                                                   num_elems, err);
      MSQ_CHKERR(err);
      
        // Get the topologies of those elements
      if (elemTopologiesSize < num_elems)
      {
        delete [] elemTopologies;
        elemTopologies = new EntityTopology[num_elems];
        (*currentMesh)->elements_get_topologies(elemArray,
                                                elemTopologies,
                                                num_elems, err);
        MSQ_CHKERR(err);
        elemTopologiesSize = num_elems;
      }
      
        // Figure out how many vertices we need to allocate
      size_t num_vert_uses = 1;
      size_t i;
      for (i = 0; i < num_elems; ++i)
      {
        num_vert_uses += vertices_in_topology(elemTopologies[i]);
      }
        // All elems share at least 1 vertex (the center vertex)
      size_t num_verts = num_vert_uses - num_elems;
      
        // Get the vertices attached to those elements
      if (vertArraySize < num_verts)
      {
        delete [] vertArray;
        vertArray = new Mesh::VertexHandle[num_verts];
        delete [] vertexOnBoundary;
        vertexOnBoundary = new bool[num_verts];
        vertArraySize = num_verts;
      }
      if (csrDataSize < num_vert_uses)
      {
        delete [] csrData;
        csrData = new size_t[num_vert_uses];
        csrDataSize = num_vert_uses;
      }
      if (csrOffsetsSize < num_elems + 1)
      {
        delete [] csrOffsets;
        csrOffsets = new size_t[num_elems + 1];
        csrOffsetsSize = num_elems + 1;
      }
      (*currentMesh)->elements_get_attached_vertices(elemArray,
                                                     num_elems,
                                                     vertArray,
                                                     num_verts,
                                                     csrData,
                                                     num_vert_uses,
                                                     csrOffsets,
                                                     err); MSQ_CHKERR(err);
      
        // Allocate the space for the vertices in the PatchData
      pd.reserve_vertex_capacity(num_verts, err); MSQ_CHKERR(err);
      
        // Get the coordinates of the vertices and its flags.
      MsqVertex* pd_vert_array = pd.get_vertex_array(err);
        //get the coordinates
      (*currentMesh)->vertices_get_coordinates(vertArray,
                                               pd_vert_array,
                                               num_verts,
                                               err); MSQ_CHKERR(err);
      for (i = 0; i < num_verts; i++)
      {
        
          // If it's not the center vertex, mark it as hard fixed
        if (vertArray[i] != vertex)
        {
            // Get its flags
          (*currentMesh)->vertex_get_byte(vertArray[i],
                                          &(pd_vert_array[i].vertexBitFlags),
                                          err); MSQ_CHKERR(err);
          pd_vert_array[i].vertexBitFlags |= MsqVertex::MSQ_HARD_FIXED;
        }
          //else it is the center vertex.  We therefore already have
          //the fixed flag stored center_fixed_byte.  The hard fixed
          //flag has already been removed (when flag was retreived).
        else{
          pd_vert_array[i].vertexBitFlags = (center_fixed_byte);
        }
        
          // Add its handle to the patch
        pd.vertexHandlesArray[i] = vertArray[i];
      }
      pd.numVertices = num_verts;
      
        // Allocate space for the elements in the PatchData
      pd.reserve_element_capacity(num_elems, err); MSQ_CHKERR(err);
      
        // Put the elements into the PatchData
      MsqMeshEntity* pd_elem_array = pd.get_element_array(err);
      for (i = 0; i < num_elems; ++i)
      {
        pd_elem_array[i].set_element_type(elemTopologies[i]);
        for (size_t j = vertices_in_topology(elemTopologies[i]); j--; )
        {
          pd_elem_array[i].set_vertex_index(j, csrData[csrOffsets[i]+j]);
        }
          // Copy the element's handle to the patch
        pd.elementHandlesArray[i] = elemArray[i];
      }
      pd.numElements = num_elems;
    }
    FUNCTION_TIMER_END();  
    return true;
    case PatchData::GLOBAL_PATCH:
    {
        // We only support global patches for a single Mesh
      if (meshSet.size() != 1)
      {
        err.set_msg("Global patches only supported for single-Mesh MeshSets.");
        FUNCTION_TIMER_END();  
        return false;
      }
      
      pd.mType = PatchData::GLOBAL_PATCH;
      
        // for a global patch, we always reset to start of the mesh.
      reset(err);
      
      size_t i;
      
        // Get all vertices
      size_t num_verts = (*currentMesh)->get_total_vertex_count(err);
      MSQ_CHKERR(err);
      if (vertArraySize < num_verts)
      {
        delete [] vertArray;
        vertArray = new Mesh::VertexHandle[num_verts];
        delete [] vertexOnBoundary;
        vertexOnBoundary = new bool[num_verts];
        vertArraySize = num_verts;
      }
      (*currentMesh)->get_all_vertices(vertArray, num_verts,
                                       err); MSQ_CHKERR(err);
      
        // Put them into the patch
      pd.reserve_vertex_capacity(num_verts, err); MSQ_CHKERR(err);
      MsqVertex* pd_vert_array = pd.get_vertex_array(err);

      (*currentMesh)->vertices_get_coordinates(vertArray,
                                               pd_vert_array,
                                               num_verts,
                                               err); MSQ_CHKERR(err);
      (*currentMesh)->vertices_are_on_boundary(vertArray, vertexOnBoundary,
                                               num_verts, err); MSQ_CHKERR(err);
      for (i = 0; i < num_verts; i++)
      {
          // Get its flags
        /*(*currentMesh)->vertex_get_byte(vertArray[i],
                                        &(pd_vert_array[i].vertexBitFlags),
                                        err); MSQ_CHKERR(err);*/
          // Set its hard-fixed flag
        if (/*(*currentMesh)->vertex_is_fixed(vertArray[i], err) ||*/
           vertexOnBoundary[i])
        {
          pd_vert_array[i].vertexBitFlags |= MsqVertex::MSQ_HARD_FIXED;
        }
        else
        {
          pd_vert_array[i].vertexBitFlags &= ~(MsqVertex::MSQ_HARD_FIXED);
        }
        MSQ_CHKERR(err);
          // Add its handle to the patch data
        pd.vertexHandlesArray[i] = vertArray[i];
      }
      pd.numVertices = num_verts;
      
        // Get all elements
      size_t num_elems = (*currentMesh)->get_total_element_count(err);
      MSQ_CHKERR(err);
      if (elemArraySize < num_elems)
      {
        delete [] elemArray;
        elemArray = new Mesh::ElementHandle[num_elems];
        elemArraySize = num_elems;
      }
      (*currentMesh)->get_all_elements(elemArray, num_elems, err); MSQ_CHKERR(err);

        // Get the topologies of those elements
      if (elemTopologiesSize < num_elems)
      {
        delete [] elemTopologies;
        elemTopologies = new EntityTopology[num_elems];
        elemTopologiesSize = num_elems;
      }
      (*currentMesh)->elements_get_topologies(elemArray, elemTopologies,
                                              num_elems,err);
     
      size_t num_attached_vtx=0;
      for (i = 0; i < num_elems; ++i)
        num_attached_vtx += vertices_in_topology(elemTopologies[i]);
      
        // Put them into the patch
      pd.reserve_element_capacity(num_elems, err); MSQ_CHKERR(err);
      MsqMeshEntity* pd_elem_array = pd.get_element_array(err);
      for (i = 0; i < num_elems; ++i)
      {
        pd_elem_array[i].set_element_type(elemTopologies[i]);
        (*currentMesh)->element_get_attached_vertex_indices(
          elemArray[i],
          pd_elem_array[i].get_modifiable_vertex_index_array(),
          MSQ_MAX_NUM_VERT_PER_ENT,
          err); MSQ_CHKERR(err);
        pd.elementHandlesArray[i] = elemArray[i];
      }
      pd.numElements = num_elems;
    }
    FUNCTION_TIMER_END();
//    pd.print(); //dbg
    return true;
    default:
      err.set_msg("no implementation for specified patch type.");
      FUNCTION_TIMER_END();  
      return false;
  }

  FUNCTION_TIMER_END();
  return true;
}

// Currently, the only thing supported is updating each vertices
// coordinates and flags.  Connectivity changes aren't supported yet.
void Mesquite::MeshSet::update_mesh(const PatchData &pd, MsqError &err)
{
  if (pd.numVertices == 0)
    return;
  
  size_t i;
  
  switch (pd.type())
  {
    // If the patch type is marked as local,
    // all handles belong to the currentMesh.
  case PatchData::ELEMENTS_ON_VERTEX_PATCH:
      // For each vertex, update the coordinates
        // and the "mesquite byte".
      for (i = 0; i < pd.numVertices; i++)
      {
        (*currentMesh)->vertex_set_coordinates(pd.vertexHandlesArray[i],
                                               pd.vertexArray[i],
                                               err); MSQ_CHKERR(err);
        (*currentMesh)->vertex_set_byte(pd.vertexHandlesArray[i],
                                        pd.vertexArray[i].vertexBitFlags,
                                        err); MSQ_CHKERR(err);
      }
      break;
      
    // If the patch type is marked as global,
    // the handles may belong to more than
    // one Mesh.
  case PatchData::GLOBAL_PATCH:
    {
      list<Mesquite::Mesh*>::iterator mesh_itr = meshSet.begin();
      Mesquite::Mesh* cur_mesh = *mesh_itr;
      Mesquite::VertexIterator *vert_itr = cur_mesh->vertex_iterator(err);
      for (i = 0; i < pd.numVertices; i++)
      {
        if (vert_itr->is_at_end())
        {
          mesh_itr++;
          cur_mesh = *mesh_itr;
          delete vert_itr;
          vert_itr = cur_mesh->vertex_iterator(err); MSQ_CHKERR(err);
        }
        cur_mesh->vertex_set_coordinates(pd.vertexHandlesArray[i],
                                         pd.vertexArray[i],
                                         err); MSQ_CHKERR(err);
        cur_mesh->vertex_set_byte(pd.vertexHandlesArray[i],
                                  pd.vertexArray[i].vertexBitFlags,
                                  err); MSQ_CHKERR(err);
      }
      delete vert_itr;
    }
    break;
  default:
    {
      err.set_msg("PatchData Type not accepted yet."); MSQ_CHKERR(err);
      break;
    }
  }
}

bool MeshSet::clear_all_soft_fixed_flags(MsqError &err)
{
    //variable to store the center vertex's fixed flag
  MsqVertex::FlagMask fixed_byte;
  bool finished_with_vertices=false;
    // initialize everything.
  if (!vertexIterator)
    reset(err);
    // currentVertex is pointing at next potential center vertex.
    
  while(!finished_with_vertices){
      // Move to next mesh if necessary
    if (vertexIterator->is_at_end())
    {
      delete vertexIterator;
      ++currentMesh;
      if (currentMesh == meshSet.end())
      {
        vertexIterator = NULL;
        finished_with_vertices=true;
      }
      if(!finished_with_vertices){
        vertexIterator = (*currentMesh)->vertex_iterator(err); MSQ_CHKERR(err);
      }
    }
      //otherwise we check to see if this vertex has been culled
    else{
        //get the fixed_bit_flag 
      (*currentMesh)->vertex_get_byte(**vertexIterator,&fixed_byte, err);
      fixed_byte &= (~MsqVertex::MSQ_SOFT_FIXED);
      (*currentMesh)->vertex_set_byte(**vertexIterator,fixed_byte, err);
      ++(*vertexIterator);
      MSQ_CHKERR(err);
    }
  }
  return true;
}
