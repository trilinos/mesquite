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


/** \file JacobianCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "JacobianCalculator.hpp"
#include "MappingFunction.hpp"
#include "MsqError.hpp"
#include "TopologyInfo.hpp"

namespace Mesquite {

void JacobianCalculator::get_derivatives( const MappingFunction* func,
                                          unsigned ho_bits,
                                          unsigned dim, unsigned num,
                                          MsqError& err )
{
  if (!func) {
    MSQ_SETERR(err)( "No mapping function.", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  unsigned edim = TopologyInfo::dimension( func->element_topology() );
  
  mIndices.clear();
  mDerivs.clear();
  switch (dim) {
    case 0:
      func->derivatives_at_corner( num, ho_bits, mIndices, mDerivs, err );
      break;
    case 1:
      func->derivatives_at_mid_edge( num, ho_bits, mIndices, mDerivs, err );
      break;
    case 2:
      if (edim != 2) {
        func->derivatives_at_mid_face( num, ho_bits, mIndices, mDerivs, err );
        break;
      }
    case 3:
      func->derivatives_at_mid_elem( ho_bits, mIndices, mDerivs, err );
      break;
  }
  MSQ_CHKERR( err );
}

void JacobianCalculator::get_Jacobian_2D( const MappingFunction* mf,
                                          unsigned ho_bits,
                                          unsigned dim, unsigned num,
                                          const Vector3D* verts,
                                          MsqMatrix<3,2>& J_out,
                                          MsqError& err )
{
  get_derivatives( mf, ho_bits, dim, num, err ); MSQ_ERRRTN(err);
  msq_std::vector<size_t>::const_iterator i = mIndices.begin();
  msq_std::vector<double>::const_iterator d = mDerivs.begin();
  Vector3D c[2] = {Vector3D(0,0,0), Vector3D(0,0,0)};
  for (; i != mIndices.end(); ++i) {
    c[0] += *d * verts[*i]; ++d;
    c[1] += *d * verts[*i]; ++d;
  }
  J_out.set_column( 0, *(MsqMatrix<3,1>*)&c[0] );
  J_out.set_column( 1, *(MsqMatrix<3,1>*)&c[1] );
}

void JacobianCalculator::get_Jacobian_3D( const MappingFunction* mf,
                                          unsigned ho_bits,
                                          unsigned dim, unsigned num,
                                          const Vector3D* verts,
                                          MsqMatrix<3,3>& J_out,
                                          MsqError& err )
{
  get_derivatives( mf, ho_bits, dim, num, err ); MSQ_ERRRTN(err);
  msq_std::vector<size_t>::const_iterator i = mIndices.begin();
  msq_std::vector<double>::const_iterator d = mDerivs.begin();
  Vector3D c[3] = {Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0)};
  for (; i != mIndices.end(); ++i) {
    c[0] += *d * verts[*i]; ++d;
    c[1] += *d * verts[*i]; ++d;
    c[2] += *d * verts[*i]; ++d;
  }
  J_out.set_column( 0, *(MsqMatrix<3,1>*)&c[0] );
  J_out.set_column( 1, *(MsqMatrix<3,1>*)&c[1] );
  J_out.set_column( 2, *(MsqMatrix<3,1>*)&c[2] );
}


} // namespace Mesquite
