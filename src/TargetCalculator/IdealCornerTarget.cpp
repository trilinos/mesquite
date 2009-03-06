/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file IdealCornerTarget.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdealCornerTarget.hpp"
#include "IdealElements.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"
#include "MsqMatrix.hpp"
#include "PatchData.hpp"
#include "ElemSampleQM.hpp"

namespace Mesquite {

bool IdealCornerTarget::surface_targets_are_3D() const
  { return true; }
 
bool IdealCornerTarget:: get_3D_target( PatchData& pd, 
                                        size_t element,
                                        const SamplePoints* pts,
                                        Sample sample,
                                        MsqMatrix<3,3>& W,
                                        MsqError& err )
{
  const MsqMeshEntity& elem = pd.element_by_index(element);
  const EntityTopology type = elem.get_element_type();

  unsigned corner = sample.number;
  if (sample.dimension != 0) {
    MSQ_SETERR(err)("IdealCornerTarget cannot generate targets at sample "
                    "locations other than corners", MsqError::INVALID_STATE);
    return false;
  }
  
  const Vector3D* coords = unit_edge_element( type, true );
  unsigned adj_count;
  const unsigned* adj = TopologyInfo::adjacent_vertices( type, corner, adj_count );
  Matrix3D& m = reinterpret_cast<Matrix3D&>(W);
  Vector3D n;
  m.set_column( 0, coords[adj[0]] - coords[corner] );
  m.set_column( 1, coords[adj[1]] - coords[corner] );
  switch (adj_count) {
    case 2: 
      n = m.column(0) * m.column(1);
      n.normalize();
      if (type == TRIANGLE)
        n *= MSQ_3RT_2_OVER_6RT_3;
      m.set_column( 2, n );
      break;
    case 3:
      m.set_column( 2, coords[adj[2]] - coords[corner] );
      break;
    default:
      MSQ_SETERR(err)("Invalid corner", MsqError::INVALID_ARG);
      return false;
  }
  return true;
}

bool IdealCornerTarget:: get_2D_target( PatchData& , 
                                        size_t ,
                                        const SamplePoints* ,
                                        Sample ,
                                        MsqMatrix<3,2>& ,
                                        MsqError& err )
{
  MSQ_SETERR(err)("IdealCornerTarget cannot be used with jacobian-based "
                  "target metrics", MsqError::INVALID_STATE);
  return false;
}

} // namespace Mesquite
