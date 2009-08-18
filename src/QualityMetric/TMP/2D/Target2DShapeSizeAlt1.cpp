/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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
 
    (2009) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file Target2DShapeSizeAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Target2DShapeSizeAlt1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

msq_std::string Target2DShapeSizeAlt1::get_name() const
  { return "ShapeSize1"; }

/** \f$ \mu(T) = \frac{|T|^2+2}{2\psi(T)} - 1 \f$
 *  \f$ \psi(T) = \sqrt{|T|^2 + 2 \tau}\f$
 *  \f$ \tau = det(T) \f$
 */
bool Target2DShapeSizeAlt1::evaluate( const MsqMatrix<2,2>& A, 
                                      const MsqMatrix<2,2>& W, 
                                      double& result, 
                                      MsqError&  )
{
  MsqMatrix<2,2> T = A * inverse(W);
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt( frob_sqr + 2.0*det(T) );

  while (!Mesquite::divide(frob_sqr+2,2*psi,result)) {
    T(0,0) += DBL_EPSILON;
    T(1,1) += DBL_EPSILON;
    frob_sqr = sqr_Frobenius(T);
    psi = sqrt( frob_sqr + 2.0*det(T) );
  }

  result -= 1.0;
  return true;
}

/*
bool Target2DShapeSizeAlt1::evaluate_with_grad( const MsqMatrix<2,2>& A, 
                                                const MsqMatrix<2,2>& W, 
                                                double& result, 
                                                MsqMatrix<2,2>& deriv_wrt_A,
                                                MsqError& err )
{
}
*/

/*
bool Target2DShapeSizeAlt1::evaluate_with_hess( const MsqMatrix<2,2>& A, 
                                                const MsqMatrix<2,2>& W, 
                                                double& result, 
                                                MsqMatrix<2,2>& deriv_wrt_A,
                                                MsqMatrix<2,2> second[3],
                                                MsqError& err )
{
}
*/

} // namespace Mesquite
