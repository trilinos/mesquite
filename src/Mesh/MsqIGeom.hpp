
/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */
/*!
  \file   MsqIGeom.hpp
  \brief  Mesquite::MeshDomain implemented on ITAPS iGeom API
  \author Jason Kraftcheck
  \date   2007-08-14
*/

#ifndef MSQ_IGEOM_HPP
#define MSQ_IGEOM_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "iGeom.h"
#include "iRel.h"

namespace MESQUITE_NS
{

    /**\brief A base class describing a Mesquite::MeshDomain implemented
     *        on top of the ITAPS iGeom and iRel APIs.
     *
     * A base class for Mesquite::MeshDomain implementations on top
     * of the ITAPS iGeom and iRel APIs and static
     * methods for constructing concrete implementations.  The concrete
     * implementations are not themselves in a header to avoid the
     * need to include all ITAPS headers if this header is included, for
     * example indirectly through Mesquite_all_headers.hpp
     */
  class MsqIGeom : public Mesquite::MeshDomain
  {
    public:
    
      /**\brief Create MeshDommain from iGeom and iRel instances
        *
        * Create an instance of an implementation of MeshDomain that uses the 
        * iRel interface to get a handle for the geometric entity 
        * associated with a mesh entity and the iGeom interface to
        * evaluate the geometry.
        */
    static MsqIGeom* create( iGeom_Instance geom, 
                             iRel_Instance irel_iface,
                             iRel_RelationHandle irel_instance,
                             MsqError& err );
    
      /**\brief Create a MeshDomain for a single geometric entity using
       *        the iGeom API for geometric evaluation.
       *
       * Create a iGeom MeshDomain for a single geometric entity.
       * This implementation will be faster than the one that uses the
       * classification interface because it assumes all entities in the
       * mesh are in the interior of the single geometric entity specified.
       * This implemenation in intended to be used only in the case where 
       * the mesh of a single surface is being smoothed and the mesh 
       * vertices on the boundary of the surface are fixed.
       */
    static MsqIGeom* create( iGeom_Instance geom,
                             iBase_EntityHandle geom_ent_handle,
                             MsqError& err );
  
    virtual ~MsqIGeom();
  };
}

#endif
