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


/** \file JacobianCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_JACOBIAN_CALCULATOR_HPP
#define MSQ_JACOBIAN_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "MsqMatrix.hpp"
#include "Vector3D.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <vector.h>
#else
# include <vector>
#endif

namespace Mesquite {

class MappingFunction;

/**\brief Calculate Jacobian matrices given vertex coordinates and MappingFunction
*/
class JacobianCalculator
{
public:    
  /**\brief Calculate Jacobian for surface element
   *
   * Calculate the Jacobian matrix at a specified location in a surface element.
   *\param mf The mapping function
   *\param ho_bits bit mask indicating which higher-order nodes are present in the element
   *               (zero for a linear element)
   *\param dim The dimension of the sub-entity the sample point lies on
   *           (0->corner, 1->mid-edge, 2->mid-face, ...)
   *\param num The number of the entity at which to sample (corner numer,
   *           edge number, etc.)
   *\param J_out The resulting Jacobian matrix.
   */
  void get_Jacobian_2D( const MappingFunction* mf,
                        unsigned ho_bits,
                        unsigned dim, unsigned num,
                        const Vector3D* vertex_coords,
                        MsqMatrix<3,2>& J_out,
                        MsqError& err );
  
  /**\brief Calculate Jacobian for volume element
   *
   * Calculate the Jacobian matrix at a specified location in a volume element.
   *\param mf The mapping function
   *\param ho_bits bit mask indicating which higher-order nodes are present in the element
   *               (zero for a linear element)
   *\param dim The dimension of the sub-entity the sample point lies on
   *           (0->corner, 1->mid-edge, 2->mid-face, ...)
   *\param num The number of the entity at which to sample (corner numer,
   *           edge number, etc.)
   *\param J_out The resulting Jacobian matrix.
   */
  void get_Jacobian_3D( const MappingFunction* mf,
                        unsigned ho_bits,
                        unsigned dim, unsigned num,
                        const Vector3D* vertex_coords,
                        MsqMatrix<3,3>& J_out,
                        MsqError& err );
private:
  void get_derivatives( const MappingFunction* mf,
                        unsigned ho_bits,
                        unsigned dim, unsigned num, 
                        MsqError& err );

  msq_std::vector<double> mDerivs;
  msq_std::vector<size_t> mIndices;
};




} // namespace Mesquite

#endif
