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


/** \file TetLagrangeShape.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TET_LAGRANGE_SHAPE_HPP
#define MSQ_TET_LAGRANGE_SHAPE_HPP

#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace Mesquite {

/**\brief Lagrange shape function for tetrahedral elements
 *
 * This class implements the MappingFunction interface, providing
 * a Lagrange shape function for triangle elements.
 *
 * \f$\vec{x}(r,s) = sum_{i=0}^{n-1} N_i(r,s) \vec{x_i}\f$
 * 
 * \f$N_1 = r(2r - 1)\f$
 * 
 * \f$N_2 = s(2s - 1)\f$
 * 
 * \f$N_3 = t(2t - 1)\f$
 *
 * \f$N_4 = u(2u - 1)\f$
 * 
 * \f$N_5 = 4rs\f$
 * 
 * \f$N_6 = 4st\f$
 * 
 * \f$N_7 = 4rt\f$
 * 
 * \f$N_8 = 4ru\f$
 * 
 * \f$N_9 = 4su\f$
 * 
 * \f$N_10 = 4tu\f$
 *
 * \f$u = 1 - r - s - t\f$
 */
class TetLagrangeShape : public MappingFunction
{
public:

  virtual 
  EntityTopology element_topology() const;

  virtual 
  void coefficients_at_corner( unsigned corner, 
                               unsigned nodebits,
                               msq_std::vector<double>& coeff_out,
                               MsqError& err ) const; 

  virtual 
  void coefficients_at_mid_edge( unsigned edge, 
                                 unsigned nodebits,
                                 msq_std::vector<double>& coeff_out,
                                 MsqError& err ) const;

  virtual 
  void coefficients_at_mid_face( unsigned face, 
                                 unsigned nodebits,
                                 msq_std::vector<double>& coeff_out,
                                 MsqError& err ) const;

  virtual 
  void coefficients_at_mid_elem( unsigned nodebits,
                                 msq_std::vector<double>& coeff_out,
                                 MsqError& err ) const;

  virtual 
  void derivatives_at_corner( unsigned corner, 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const;

  virtual 
  void derivatives_at_mid_edge( unsigned edge, 
                                unsigned nodebits,
                                msq_std::vector<size_t>& vertex_indices_out,
                                msq_std::vector<double>& d_coeff_d_xi_out,
                                MsqError& err ) const;

  virtual 
  void derivatives_at_mid_face( unsigned face, 
                                unsigned nodebits,
                                msq_std::vector<size_t>& vertex_indices_out,
                                msq_std::vector<double>& d_coeff_d_xi_out,
                                MsqError& err ) const;

  virtual 
  void derivatives_at_mid_elem( unsigned nodebits,
                                msq_std::vector<size_t>& vertex_indices_out,
                                msq_std::vector<double>& d_coeff_d_xi_out,
                                MsqError& err ) const;
};



} // namespace Mesquite

#endif
