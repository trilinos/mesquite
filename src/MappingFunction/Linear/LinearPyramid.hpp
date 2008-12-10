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

/** \file LinearPyramid.hpp
 *  \brief mapping function for linear prism
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_LINEAR_PYRAMID_HPP
#define MSQ_LINEAR_PYRAMID_HPP

#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace Mesquite {

/**\brief Linear mapping function for a pyramid element
 *
 * \f$\vec{x}(\vec{\xi})=\sum_{i=0}^5 N_i(\vec{\xi})\vec{x_i}\f$
 * - \f$N_0(\vec{\xi})=\frac{1}{8}(1-\xi)(1-\eta)(1-\zeta)\f$
 * - \f$N_1(\vec{\xi})=\frac{1}{8}(1+\xi)(1-\eta)(1-\zeta)\f$
 * - \f$N_2(\vec{\xi})=\frac{1}{8}(1+\xi)(1+\eta)(1-\zeta)\f$
 * - \f$N_3(\vec{\xi})=\frac{1}{8}(1-\xi)(1+\eta)(1-\zeta)\f$
 * - \f$N_4(\vec{\xi})=\frac{1}{2}(1+\zeta)\f$
 *
 * Where the mid-edge coordinates are:
 * - \f$\vec{\xi}_{e0}=( 0,-1,-1)\f$
 * - \f$\vec{\xi}_{e1}=( 1, 0,-1)\f$
 * - \f$\vec{\xi}_{e2}=( 0, 1,-1)\f$
 * - \f$\vec{\xi}_{e3}=(-1, 0,-1)\f$
 * - \f$\vec{\xi}_{e4}=(-1,-1, 0)\f$
 * - \f$\vec{\xi}_{e5}=( 1,-1, 0)\f$
 * - \f$\vec{\xi}_{e6}=( 1, 1, 0)\f$
 * - \f$\vec{\xi}_{e7}=(-1, 1, 0)\f$
 *
 * and mid-face the coordinates are:
 * - \f$\vec{\xi}_{f0}=( 0,-1,\frac{1}{3})\f$
 * - \f$\vec{\xi}_{f1}=( 1, 0,\frac{1}{3})\f$
 * - \f$\vec{\xi}_{f2}=( 0, 1,\frac{1}{3})\f$
 * - \f$\vec{\xi}_{f3}=(-1, 0,\frac{1}{3})\f$
 * - \f$\vec{\xi}_{f4}=( 0, 0, -1)\f$
 *
 * and the mid-element coorindates are
 * - \f$\vec{\xi}_m=( 0, 0,-\frac{1}{2})\f$
 */ 
class LinearPyramid : public MappingFunction
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
