/*!
  \file   MeshSet.hpp
  \brief  

  The MeshSet class provides the control interface for a certain mesh set.

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_MeshSet_hpp 
#define Mesquite_MeshSet_hpp

#ifdef HAVE_STANDARD_INCLUDES
#include <cstddef>
#else
#include <stddef.h>
#endif

#include <list>
#include <vector>
#include <iterator>
#include <string>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "PatchData.hpp"

#include "TSTT_Base.h"


namespace Mesquite
{
  class PatchData;

    /*! \class MeshSet
     !!! 
     This is a very preliminary version intended only
     to give a minimum functionality
     !!! */ 
  class MeshSet
  {
  public:
    
    MeshSet();
    ~MeshSet();
    void add_mesh(TSTT::Mesh_Handle mesh_handle, MsqError &err);
    void get_meshes(std::list<TSTT::Mesh_Handle> &mesh_list) const
      { mesh_list = meshSet; }
    
      // This should normally be set by the MeshSet itself when you
      // call add_mesh().
    void set_space_dim(int dim)
      { spaceDim = dim; }
    
    int space_dim() const
      { return spaceDim; }
    
    enum PatchType
    {
      UNDEFINED_PATCH_TYPE,
      VERTICES_ON_VERTEX_PATCH,
      ELEMENTS_ON_VERTEX_PATCH,
      GLOBAL_PATCH
    };
    
      // Tells the MeshSet what kind of data the patches should include,
      // what they should be based on.
      // The meaning of patch_param1 and patch_param2 depends on
      // patch_type.
    bool set_patch_type(PatchType patch_type,
                        int patch_param1=0,
                        int patch_param2=0);

    PatchType get_patch_type() {return mType;}
    
    /*! Sets the name of the tag that identifies fixed vertices
      within that particular MeshSet .
        By default, the tag name for fixed vertices is "fixed".
     */
    void set_fixed_vertex_tag(std::string tag_name)
    { fixedVertexTagName = tag_name; }
    
    void copy_culling_method_bits (long unsigned int bits)
    { cullingMethodBits = bits; }

    // Gets the next patch, according to what has been set in
      // set_patch_type().
    bool get_next_patch(PatchData &pd, MsqError &err);

    void reset(MsqError &err);

      // These should eventually be handled by get_next_patch()
    bool get_next_element_group(PatchData &pd, MsqError &err);
    bool get_next_node_group(PatchData &pd, MsqError &err);
    
    struct EntityEntry
    {
      TSTT::Mesh_Handle mesh;
      TSTT::Entity_Handle entity;
      EntityEntry( TSTT::Mesh_Handle m, TSTT::Entity_Handle e ) :
        mesh(m), entity(e)
        {}
    };
    
  private:
    bool get_next_vertices_set(MsqError &err);
    
    PatchType mType;
    int mParam1, mParam2;
    long unsigned int cullingMethodBits;
    std::string fixedVertexTagName; 
    std::list<TSTT::Mesh_Handle> meshSet;
    std::list<TSTT::Mesh_Handle>::iterator currentMesh;
    std::vector<EntityEntry> verticesSet;
    std::vector<EntityEntry>::iterator currentVertex;
    
    int spaceDim;
    //! TSTT::FACE or TSTT::REGION.
    //! Must be the same for all meshes added with add_mesh().
    enum TSTT::EntityType elementType;
  };


    // *********Call MsqMeshEntity::vertex_count() instead**********
// #undef __FUNC__
// #define __FUNC__ "number_of_vertices"
//   inline int number_of_vertices(enum TSTT::EntityTopology topo, MsqError &err)
//   {
//     switch(topo) {
//     case TSTT::TRIANGLE:
//       return 3;
//     case TSTT::QUADRILATERAL:
//       return 4;
//     case TSTT::TETRAHEDRON:
//       return 4;
//     case TSTT::HEXAHEDRON:
//       return 8;
//     case TSTT::PRISM:
//       return 6;
//     case TSTT::PYRAMID:
//       return 5;
//     case TSTT::SEPTAHEDRON:
//       return 7;
      
//     default:
//       err.set_msg("invalid element type.");
//       return 0;
//     }
//   }
  
  // -********** AOMD tmp TEST **********
  void test_aomd(void);

} //namespace

#endif
