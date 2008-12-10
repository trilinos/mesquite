/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_LINEAR_QUADRILATERAL_HPP
#define MSQ_LINEAR_QUADRILATERAL_HPP

/** \file LinearQuadrilateral.hpp
 *  \brief Linear mapping funtion for Quadrilateral elements.
 *  \author Jason Kraftcheck
 */
 
#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace Mesquite {

/**\brief Linear shape function for quadrilateral elements
 *
 * This class implements the MappingFunction interface, providing
 * a Linear shape function for quadrilateral elements.
 *
 * \f$\vec{x}(\xi,\eta) = \sum_{i=0}^{3} N_i(\xi,\eta) \vec{x_i}\f$
 * 
 * \f$N_0(\xi,\eta) = \frac{1}{4}(1-\xi)(1-\eta)\f$
 * 
 * \f$N_1(\xi,\eta) = \frac{1}{4}(1+\xi)(1-\eta)\f$
 * 
 * \f$N_2(\xi,\eta) = \frac{1}{4}(1+\xi)(1+\eta)\f$
 * 
 * \f$N_3(\xi,\eta) = \frac{1}{4}(1-\xi)(1+\eta)\f$
 */
class LinearQuadrilateral : public MappingFunction
{
public:

  virtual 
  EntityTopology element_topology() const;

  virtual 
  void coefficients_at_corner( unsigned corner, 
                               unsigned nodebits,
                               double* coeff_out,
                               size_t& num_coeff,
                               MsqError& err ) const; 

  virtual 
  void coefficients_at_mid_edge( unsigned edge, 
                                 unsigned nodebits,
                                 double* coeff_out,
                                 size_t& num_coeff,
                                 MsqError& err ) const;

  virtual 
  void coefficients_at_mid_face( unsigned face, 
                                 unsigned nodebits,
                                 double* coeff_out,
                                 size_t& num_coeff,
                                 MsqError& err ) const;

  virtual 
  void coefficients_at_mid_elem( unsigned nodebits,
                                 double* coeff_out,
                                 size_t& num_coeff,
                                 MsqError& err ) const;

  virtual 
  void derivatives_at_corner( unsigned corner, 
                              unsigned nodebits,
                              size_t* vertex_indices_out,
                              double* d_coeff_d_xi_out,
                              size_t& num_vtx,
                              MsqError& err ) const;

  virtual 
  void derivatives_at_mid_edge( unsigned edge, 
                                unsigned nodebits,
                                size_t* vertex_indices_out,
                                double* d_coeff_d_xi_out,
                                size_t& num_vtx,
                                MsqError& err ) const;

  virtual 
  void derivatives_at_mid_face( unsigned face, 
                                unsigned nodebits,
                                size_t* vertex_indices_out,
                                double* d_coeff_d_xi_out,
                                size_t& num_vtx,
                                MsqError& err ) const;

  virtual 
  void derivatives_at_mid_elem( unsigned nodebits,
                                size_t* vertex_indices_out,
                                double* d_coeff_d_xi_out,
                                size_t& num_vtx,
                                MsqError& err ) const;
};

} // namespace Mesquite

#endif
