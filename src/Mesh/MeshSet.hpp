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
  \file   MeshSet.hpp
  \brief  

  The MeshSet class provides the control interface for a certain mesh set.

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_MeshSet_hpp 
#define Mesquite_MeshSet_hpp

#ifdef USE_C_PREFIX_INCLUDES
#include <cstddef>
#include <cstdio>
#else
#include <stddef.h>
#include <stdio.h>
#endif

#include <list>
#include <vector>
#include <iterator>
#include <string>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "PatchDataUser.hpp"
#include "MeshInterface.hpp"

MSQ_USE(list);


namespace Mesquite
{
  class PatchData;
  class PatchDataUser;
  class PatchDataParameters;
  
  /*! \class MeshSet
    
      \brief The MeshSet class stores one or more Mesquite::Mesh pointers
       and manages access to the mesh information.
      
       MeshSet objects are passed to the various Mesquite algorithms in order
       to assess the quality, improve the mesh, etc... 
  */ 
  class MeshSet
  {
  public:
    
    MeshSet();
    ~MeshSet();
    
      //! adds a mesh to the MeshSet. 
      /*! If several meshes are added, the mesh information will be retrieved
        seamlessly, as if dealing with a single mesh. */
    void add_mesh(Mesquite::Mesh* mesh, MsqError &err);


      //! Sets the geometrical domain for the MeshSet. This can
      //! only be used with surface meshes. 
    void set_domain_constraint(MeshDomain* domain, MsqError &err);

    //! Returns the domain associated with the MeshSet from which the Patch originates.
    //! NULL if no domain is defined.
    Mesquite::MeshDomain* get_domain_constraint()
      { return mDomain; }
    
      //! returns the list of mesh pointers previously added. 
    void get_meshes(list<Mesquite::Mesh*> &mesh_list) const
      { mesh_list = meshSet; }

      //! Returns the number of coordinates in the Mesh's
      //! geometric coordinate system.
      // We want to rename this to geometric_dimension() const
    int space_dim() const
      { return spaceDim; }
    
      //! Gets the next PatchData.
      /*! The type of the patch is usually set on the algorithm with
        set_patch_type() and propagated to the MeshSet.
        This version of the get_next_patch() function is the most often used.
        It actually delegates to the original get_next_patch function, which
        has a slightly different signature. */
    bool get_next_patch(PatchData &pd,
                        PatchDataUser* pd_user,
                        MsqError &err)
      {
        return get_next_patch(pd, pd_user->get_all_parameters(), err);
      }
    
      /*! This version of get_next_patch() is rarely used, but this is
        where the implementation actually is.
        See the most frequently used signature:
        get_next_patch(PatchData &pd, PatchDataUser* pd_user, MsqError &err).
      */
    bool get_next_patch(PatchData &pd,
                        PatchDataParameters &pd_params,
                        MsqError &err);
    
      //! Resets the MeshSet object so that get_next_patch() will restart
      //! its iterations at the first vertex.
    void reset(MsqError &err);

      //! Updates the coordinates in the underlying mesh
      //! with the coordinates stored in PatchData
    void update_mesh(const PatchData &pd, MsqError &err);

      //! Sets the cullFlag.  This flag is used to dermine
      //! which vertices should be culled.  For local schemes,
      //! get_next_patch() will not build a patch around a
      //! a culled vertex.
    void set_flag_to_cull(MsqVertex::FlagMask f_m)
      {
        cullFlag=f_m;
      }
    
    bool clear_all_soft_fixed_flags(MsqError &err);

    void write_vtk(const char* out_filebase, MsqError &err);
    
    void write_gnuplot(const char* out_filebase, MsqError &err);
    
  private:
    
      //! Meshes in this MeshSet
    list<Mesquite::Mesh*> meshSet;
      //! Keeps track of which Mesh* we're currently
      //! working with in get_next_patch().
    list<Mesquite::Mesh*>::iterator currentMesh;
      //! Keeps track of where we are in the current mesh's vertex list
    Mesquite::VertexIterator *vertexIterator;
      //! The number of coordinates in this mesh (2D or 3D)
    int spaceDim;
    
      //! The topological dimension of the elements in this MeshSet,
      //! where Mesquite::TRIANGLE indicates 2D elements (faces) and
      //! Mesquite::TETRAHEDRON indicates 3D elements (regions).
      //! Must be the same for all meshes added with add_mesh().
    Mesquite::EntityTopology elementType;

      //! These are arrays that we cache so we don't have to reallocate
      //! at every patch.
    size_t *csrOffsets;
    size_t *csrData;
    Mesh::VertexHandle *vertArray;
    Mesh::ElementHandle *elemArray;
    EntityTopology *elemTopologies;
    bool *vertexOnBoundary;
    size_t csrOffsetsSize;
    size_t csrDataSize;
    size_t vertArraySize;
    size_t elemArraySize;
    size_t elemTopologiesSize;

      // This is the domain we snap everything back to.
    Mesquite::MeshDomain *mDomain;

      //Flag to tell the MeshSet which MSQ culling flag to use
      // to determine whether to build a local patch around a vertex
    Mesquite::MsqVertex::FlagMask cullFlag;

  };
  
    // -********** AOMD tmp TEST **********
  void test_aomd(void);
  
} //namespace

#endif
