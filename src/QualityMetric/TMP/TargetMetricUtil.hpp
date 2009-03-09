/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TargetMetricUtil.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_METRIC_UTIL_HPP
#define MSQ_TARGET_METRIC_UTIL_HPP

#include "Mesquite.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
#endif

namespace Mesquite {

template <unsigned R, unsigned C> class MsqMatrix;
class PatchData;
class MsqError;

/**\brief Calculate R and Z such that \f$W\prime = Z^{-1} W\f$ and 
 *        \f$A\prime = (RZ)^{-1} A\f$
 *
 * Calculate the matrices required to transform the active and target
 * matrices from the 3x2 surface domain to a 2x2 2D domain.
 *\param A    Input: Element Jacobian matrix.
 *\param W_32 Input: Target Jacobian matrix.
 *\param W_22 Output: 2D Target matrix.
 *\param RZ   Output: Product of R and Z needed to calculate the 2D 
 *            element matrix.
 */
void surface_to_2d( const MsqMatrix<3,2>& A, 
                    const MsqMatrix<3,2>& W_32,
                    MsqMatrix<2,2>& W_22,
                    MsqMatrix<3,2>& RZ );
/*
void surface_to_2d( const MsqMatrix<3,2>& A_in,
                    const MsqMatrix<3,2>& W_in,
                    MsqMatrix<2,2>& A_out,
                    MsqMatrix<2,2>& W_out );
*/
void get_sample_pt_evaluations( PatchData& pd,
                                msq_std::vector<size_t>& handles,
                                bool free,
                                MsqError& err );
                    
void get_elem_sample_points( PatchData& pd,
                             size_t elem,
                             msq_std::vector<size_t>& handles,
                             MsqError& err );
                    
} // namespace Mesquite

#endif
