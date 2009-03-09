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


/** \file IdentityTarget.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdentityTarget.hpp"
#include "MsqMatrix.hpp"

namespace Mesquite {
 
bool IdentityTarget::get_3D_target( PatchData& , 
                                    size_t ,
                                    const SamplePoints* ,
                                    Sample ,
                                    MsqMatrix<3,3>& W,
                                    MsqError&  )
{
  W = MsqMatrix<3,3>(1.0);
  return true;
}

bool IdentityTarget::get_2D_target( PatchData& , 
                                    size_t ,
                                    const SamplePoints* ,
                                    Sample ,
                                    MsqMatrix<3,2>& W,
                                    MsqError&  )
{
  W = MsqMatrix<3,2>(1.0);  
  return true;
}

} // namespace Mesquite
