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


/** \file DefaultCornerTarget.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "DefaultCornerTarget.hpp"
#include "MsqError.hpp"
#include "MsqMatrix.hpp"
#include "PatchData.hpp"
#include "ElemSampleQM.hpp"

namespace Mesquite {

bool DefaultCornerTarget::surface_targets_are_3D() const
  { return true; }
 
bool DefaultCornerTarget:: get_3D_target( PatchData& pd, 
                                        size_t element,
                                        const SamplePoints* pts,
                                        Sample sample,
                                        MsqMatrix<3,3>& W,
                                        MsqError& err )
{
  const EntityTopology type = pd.element_by_index(element).get_element_type();
  static const double SIXTH_ROOT_OF_TWO = msq_stdc::pow(2., 1./6.);

  if (sample.dimension != 0) {
    MSQ_SETERR(err)("IdealCornerTarget cannot generate targets at sample "
                    "locations other than corners", MsqError::INVALID_STATE);
    return false;
  }

  switch( type )
  {
    case TRIANGLE:
    {
#ifndef MSQ_JANUS
      const double v_tri[] = {1., 0.5, 0., 0., MSQ_SQRT_THREE/2., 0., 0., 0., 1.};
#else    
      const double v_tri[] = {1., 0.5, 0., 0., 0.8660254037844385965883021 , 0., 0., 0., 1.};
#endif    
      W.set(v_tri);
      return true;
    }

    case QUADRILATERAL:
    case HEXAHEDRON:
    {
      const double ident[] = { 1.0, 0.0, 0.0, 
                               0.0, 1.0, 0.0, 
                               0.0, 0.0, 1.0 };
      W.set(ident);
      return true;
    }

    case TETRAHEDRON:
    {
#ifndef MSQ_JANUS
      const double v_tet[] = {1., 0.5, 0.5, 0., MSQ_SQRT_THREE/2., MSQ_SQRT_THREE/6., 0., 0., MSQ_SQRT_TWO/MSQ_SQRT_THREE};
#else 
      const double v_tet[] = {1., 0.5, 0.5, 0.,  0.8660254037844385965883021, 0.2886751345948128655294340, 0., 0., 0.8164965809277261454823815};
#endif    
      W.set(v_tet);
      W *= SIXTH_ROOT_OF_TWO;
      return true;
    }

    case PYRAMID:
    {
      const double v_pyr[] = { 1.0, 0.0, 0.5, 
                               0.0, 1.0, 0.5, 
                               0.0, 0.0, 1.0 };
      W.set(v_pyr);
      return true;
    }

    case PRISM:
    {
      const double s = MSQ_3RT_2_OVER_6RT_3;
      const double h = MSQ_3RT_2_OVER_6RT_3 * MSQ_SQRT_THREE_DIV_TWO;
      const double v_wdg[] = {  s, 0.5*s, 0,  
                                0,  h,    0,
                                0,  0,    s };
      W.set(v_wdg);
      return true;
    } 

    default:
      MSQ_SETERR(err)( MsqError::UNSUPPORTED_ELEMENT,
                "Unsupported element type (%d)\n",
                (int)type);
      return false;
  }
}

bool DefaultCornerTarget:: get_2D_target( PatchData& , 
                                        size_t ,
                                        const SamplePoints* ,
                                        Sample ,
                                        MsqMatrix<3,2>& ,
                                        MsqError& err )
{
  MSQ_SETERR(err)("DefaultCornerTarget cannot be used with jacobian-based "
                  "target metrics", MsqError::INVALID_STATE);
  return false;
}

} // namespace Mesquite
