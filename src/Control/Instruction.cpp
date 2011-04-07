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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file Instruction.cpp
Default implementation for Instruction::loop_over_mesh( ParallelMesh* mesh, ...)
so we do not force every class that derives from Instruction to implement this
parallel function.

  \author Martin Isenburg
  \date   2008-04-03
 */

#include "Instruction.hpp"
#include "MeshInterface.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"
#include "Settings.hpp"

using namespace Mesquite;

Instruction::~Instruction()
{}

double Instruction::loop_over_mesh( ParallelMesh* mesh, 
				    MeshDomain* domain, 
				    const Settings* settings,
				    MsqError& err )
{
  return loop_over_mesh((Mesh*)mesh, domain, settings, err);
}

void Instruction::initialize_vertex_byte( Mesh* mesh,
                                          MeshDomain* domain,
                                          const Settings* settings,
                                          MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err ); MSQ_ERRRTN(err);
  std::vector<unsigned char> bytes( verts.size(), 0 );
  std::vector<bool> flags;
  
    // Normally we start out with all bits cleared.  However, if
    // we're doing SLAVE_CALCULATED mode then we need to perserve
    // the exinsting value of that bit because it was presumably
    // set by a previous tool in the InstructionQueue.
  if (settings && settings->get_slaved_ho_node_mode() == Settings::SLAVE_CALCULATED) {
    mesh->vertices_get_byte( arrptr(verts), arrptr(bytes), verts.size(), err ); MSQ_ERRRTN(err);
    for (size_t i = 0; i < bytes.size(); ++i)
      bytes[i] &= (unsigned char)MsqVertex::MSQ_DEPENDENT;
  }
    // If getting slaved flag for higher-order nodes from application,
    // copy it into vertex byte now.
  else if (settings && settings->get_slaved_ho_node_mode() == Settings::SLAVE_FLAG) {
    mesh->vertices_get_slaved_flag( arrptr(verts), flags, verts.size(), err ); MSQ_ERRRTN(err);
    for (size_t i = 0; i < bytes.size(); ++i)
      if (flags[i])
        bytes[i] |= (unsigned char)MsqVertex::MSQ_DEPENDENT;
  }
    // If slave mode is SLAVE_ALL, that needs to be handled later
    // because we need to check things on a per-element basis.
    // For now, leave it to be handled in PatchData because we
    // already have the element array there.
  
    // If using application-specified fixed flag, copy that into
    // vertex byte now.
  if (!settings || settings->get_fixed_vertex_mode() == Settings::FIXED_FLAG) {
    mesh->vertices_get_fixed_flag( arrptr(verts), flags, verts.size(), err ); MSQ_ERRRTN(err);
    for (size_t i = 0; i < bytes.size(); ++i)
      if (flags[i])
        bytes[i] |= (unsigned char)MsqVertex::MSQ_HARD_FIXED;
  }
    // otherwise need to mark vertices as fixed based on geometric domain
  else {
    if (!domain) {
      MSQ_SETERR(err)("Request to fix vertices by domain classification "
                      "requres a domain.", MsqError::INVALID_STATE );
      return;
    }
    const unsigned short dim = settings->get_fixed_vertex_mode();
    assert(dim < 4u);
    std::vector<unsigned short> dof( verts.size() );
    domain->domain_DoF( arrptr(verts), arrptr(dof), verts.size(), err ); MSQ_ERRRTN(err);
    for (size_t i = 0; i < bytes.size(); ++i)
      if (dof[i] <= dim)
        bytes[i] |= (unsigned char)MsqVertex::MSQ_HARD_FIXED;
  }

    // Copy flag values back to mesh instance
  mesh->vertices_set_byte( arrptr(verts), arrptr(bytes), verts.size(), err ); MSQ_ERRRTN(err);
}

  
