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
  \file   MeshTSTT.hpp
  \brief  Provide interface to TSTT mesh database
  
          Provide interface to TSTT mesh database and iterators
          for a subset of the database's  mesh.


  \author Jason Kraftcheck
  \date   2004-10-22
*/

#ifndef MESQUITE_MESH_TSTT_HPP
#define MESQUITE_MESH_TSTT_HPP

#include "MeshInterface.hpp"

namespace TSTTM {
  class Mesh;
}

namespace Mesquite
{
  /*!  \class MeshTSTT

  \brief  Provide interface to TSTT mesh database
  
          Provide interface to TSTT mesh database and iterators
          for a subset of the database's  mesh.


  \author Jason Kraftcheck
  \date   2004-10-22
*/
  class MeshTSTT : public Mesquite::Mesh
  {
  public:

    /** \brief factory method 
     *
     * Create instance of MeshTSTT.  
     * If \ref set_active_set is not called on resulting object,
     * will default to the ENTIRE mesh.
     */
    static MeshTSTT* create( TSTTM::Mesh& mesh, void* meshset, MsqError& err );
    
    virtual ~MeshTSTT();
    
      /** \brief set mesh to be smoothed.
       *
       * Set the mesh which Mesquite is to smooth.  Optionally
       * specify fixed vertices.
       * NOTE: If an active set is not specified, the default
       *       is to use the global set (the ENTIRE mesh.)
       *
       *\param element_set TSTT entity set handle for set containing
       *                  mesh elements for which quality is to be improved.
       *                  Any non-element entities in the set will be
       *                  ignored.
       *\param fixed_vertices TSTT entity set handle containing list of 
       *                  vertices which are not to be moved.  If NULL,
       *                  it is assumed that all vertices may be moved.
       */
    virtual void set_active_set( void* element_set,  void* fixed_vertices, MsqError&  ) = 0;
  };
}

#endif // MESQUITE_MESH_TSTT_HPP
