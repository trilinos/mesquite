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

#ifndef MSQ_MAPPING_FUNCTION_HPP
#define MSQ_MAPPING_FUNCTION_HPP

/** \file MappingFunction.hpp
 *  \brief Header containg defintion of MappingFunction
 *  \author Jason Kraftcheck
 */


#include "Mesquite.hpp"
#include <vector>

namespace Mesquite {

class MsqError;

/**\brief An interface for a mapping function of the form 
 * \f$\vec{x}(\vec{\xi})=\sum_{i=1}^n N_i(\vec{\xi})\vec{x_i}\f$,
 * where \f$\vec{x_i}\f$ is a point
 * in \f$\mathbf{R}^3\f$ (i.e. \f$x_i,y_i,z_i\f$),
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \end{array}\right\}\f$ 
 * for surface elements and 
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \zeta_i\\ \end{array}\right\}\f$ 
 * for volume elements. 
 *
 * This is an interface for describing a mapping function for a
 * single element topology.  A mapping function is assumed to be
 * of the following form: 
 * \f$\vec{x}(\vec{\xi})=\sum_{i=1}^n N_i(\vec{\xi})\vec{x_i}\f$
 * where \f$n\f$ is the number of nodes in the element,
 * \f$\vec{x_i}\f$ is a point
 * in \f$\mathbf{R}^3\f$ (i.e. \f$x_i,y_i,z_i\f$), and
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \end{array}\right\}\f$ 
 * for surface elements and 
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \zeta_i\\ \end{array}\right\}\f$ 
 * for volume elements.  For example,
 * for a linear quadrilateral element, the mapping function will be
 * of the form: 
 * \f$\vec{x}(\xi,\eta)=N_1(\xi,\eta)\vec{x_1}
 *                     +N_2(\xi,\eta)\vec{x_2}
 *                     +N_3(\xi,\eta)\vec{x_3}
 *                     +N_4(\xi,\eta)\vec{x_4}\f$
 *
 * A single implementation of this interface may support multiple 
 * types of elements of the same topology.  Element types within
 * a topology may vary by the presences or lack there of of mid-edge,
 * mid-face, and mid-element nodes.
 */
class MappingFunction 
{
public:

  virtual
  ~MappingFunction() {}

  /**\brief Get Mesquite::EntityTopology handled by this mapping function */
  virtual 
  EntityTopology element_topology() const = 0;

  /**\name Mapping Function Coefficients
   * This group of methods return the list of scalar values (\f$N_i\f$'s) resulting
   * from the evaluation of the mapping function coefficient terms
   * \f$N_1(\vec{\xi}), N_2(\vec{\xi}), \ldots, N_n(\vec{\xi})\f$
   * for a given \f$\vec{\xi}\f$.
   *\param nodebits This is a list of which mid-nodes are present in the 
   *                element.  The sequence is the same as the canonical 
   *                ordering of the mid-nodes in the form of the element with
   *                all possible mid-nodes, beginning with the least signficiant
   *                bit.  A 1-bit means that the corresponding
   *                higher-order node is present in the element.  This value
   *                of nodebits is zero for a linear element.  
   *                For polygons, the number of nodes in the polygon should
   *                be passed in this argument.  Polygons are not allowed
   *                higher-order nodes.
   *\param coefficients_out The coefficients (\f$N_i(\vec{\xi})\f$) for each 
   *                vertex in the element, except for the higher-order nodes 
   *                for which the corresponding bit in nodebits is zero.
   */
  /*@{*/

  /**Get mapping function coefficients for an \f$\vec{\xi}\f$ corresponding
   * to the corner of the element.
   *\param corner The element corner, specified as the position of the
   *              corresponding vertex in the canonical ordering of the
   *              element.
   */
  virtual 
  void coefficients_at_corner( unsigned corner, 
                               unsigned nodebits,
                               msq_std::vector<double>& coeff_out,
                               MsqError& err ) const = 0; 

  /**Get mapping function coefficients for an \f$\vec{\xi}\f$ corresponding
   * to the middle of an edge of the element.
   *\param edge The element edge, specified using the canoncial ordering
   *            of edges for the element topology.
   */
  virtual 
  void coefficients_at_mid_edge( unsigned edge, 
                                 unsigned nodebits,
                                 msq_std::vector<double>& coeff_out,
                                 MsqError& err ) const = 0;

  /**Get mapping function coefficients for an \f$\vec{\xi}\f$ corresponding
   * to the middle of an face of the element.
   *\param face The element face, specified using the canoncial ordering
   *            of faces for the element topology.
   */
  virtual 
  void coefficients_at_mid_face( unsigned face, 
                                 unsigned nodebits,
                                 msq_std::vector<double>& coeff_out,
                                 MsqError& err ) const = 0;

  /**Get mapping function coefficients for an \f$\vec{\xi}\f$ corresponding
   * to the center of the element.
   */
  virtual 
  void coefficients_at_mid_elem( unsigned nodebits,
                                 msq_std::vector<double>& coeff_out,
                                 MsqError& err ) const = 0;
  /*@}*/

  /**\name Mapping Function Derivatives
   * This group of methods return the partial derivatives of the mapping
   * function coefficient terms
   * \f$\nabla N_1(\vec{\xi}), \nabla N_2(\vec{\xi}), \ldots, \nabla N_n(\vec{\xi})\f$
   * evaluated for a given \f$\vec{\xi}\f$, where \f$\vec{x_i}\f$ is a point
   * in \f$\mathbf{R}^3\f$ (i.e. \f$x_i,y_i,z_i\f$).  
   * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \end{array}\right\}\f$ 
   * for surface elements and 
   * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \zeta_i\\ \end{array}\right\}\f$ 
   * for volume elements.
   *
   * The list of returned partial derivatives may be considered list of elements 
   * of a matrix \f$\mathbf{D}\f$ in row major order.  For surface elements,
   * \f$\mathbf{D}\f$ is a \f$n\times 2\f$ matrix and for volume elements it
   * is a \f$n \times 3\f$ matrix.  Each row of 
   * \f$\mathbf{D}\f$ corresponds to one of the
   * coefficient functions \f$N_i(\vec{\xi})\f$ and each column corresponds
   * to one of the components of \f$\vec{\xi}\f$ 
   * that the corresponding coefficient function is differentiated with
   * respect to. 
   *
   * \f$ \mathbf{D} = \left[ \begin{array}{ccc}
   *     \frac{\delta N_1}{\delta \xi} & \frac{\delta N_1}{\delta \eta} & \ldots \\
   *     \frac{\delta N_2}{\delta \xi} & \frac{\delta N_2}{\delta \eta} & \ldots \\
   *     \vdots & \vdots & \ddots \end{array} \right]\f$
   *
   * The Jacobian matrix (\f$\mathbf{J}\f$) of the mapping function can be calculated
   * as follows. Define a matrix \f$\mathbf{X}\f$ such that each column contains
   * the coordinates of the element nodes.
   *
   * \f$ \mathbf{X} = \left[ \begin{array}{ccc}
   *                   x_1 & x_2 & \ldots \\
   *                   y_1 & y_2 & \ldots \\
   *                   z_1 & z_2 & \ldots 
   *                  \end{array}\right]\f$
   *
   * The Jacobian matrix is then:
   *
   * \f$\mathbf{J} = \mathbf{X} \times \mathbf{D}\f$
   *
   * \f$\mathbf{X}\f$ is always \f$3\times n\f$, so \f$\mathbf{J}\f$ is
   * either \f$3\times 2\f$ (surface elements) or \f$3\times 3\f$ (volume
   * elements) depending on the dimensions of \f$\mathbf{D}\f$.
   *
   * If the Jacobian matrix of the mapping function is considered as a 
   * function of the element vertex coordinates \f$\mathbf{J}(\vec{x_1},\vec{x_2},\ldots)\f$ 
   * with \f$\vec{\xi}\f$ constant, then the gradient of that Jacobian matrix 
   * function (with respect
   * to the vertex coordinates) can be obtained from the same output list of
   * partial deravitves.
   *
   * \f$\frac{\delta \mathbf{J}}{\delta x_i} = 
   *         \left[ \begin{array}{ccc}
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots \\
   *         0 & 0 & \ldots \\ 
   *         0 & 0 & \ldots 
   *         \end{array} \right]\f$
   * \f$\frac{\delta \mathbf{J}}{\delta y_i} = 
   *         \left[ \begin{array}{ccc}
   *         0 & 0 & \ldots \\ 
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots \\
   *         0 & 0 & \ldots 
   *         \end{array} \right]\f$
   * \f$\frac{\delta \mathbf{J}}{\delta z_i} = 
   *         \left[ \begin{array}{ccc}
   *         0 & 0 & \ldots \\ 
   *         0 & 0 & \ldots \\
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots 
   *         \end{array} \right]\f$
   * 
   *
   *\param nodebits This is a list of which mid-nodes are present in the 
   *                element.  The sequence is the same as the canonical 
   *                ordering of the mid-nodes in the form of the element with
   *                all possible mid-nodes, beginning with the least signficiant
   *                bit.  A 1-bit means that the corresponding
   *                higher-order node is present in the element.  This value
   *                of nodebits is zero for a linear element.  
   *                For polygons, the number of nodes in the polygon should
   *                be passed in this argument.  Polygons are not allowed
   *                higher-order nodes.
   *\param vertices_out The list of vertices for which the corresponding
   *                coefficient in the mapping function is non-zero.  The
   *                vertices are specified by their index in the canonical
   *                ordering for an element with all mid-nodes present (i.e.
   *                first all the corner nodes, then the mid-edge nodes, ...).
   *\param d_coeff_d_xi_out The mapping function is composed of a series of 
   *                coefficient functions \f$N_i(\vec{\xi})\f$, one correspoding
   *                to the position \f$\vec{x_i}\f$ of each node in the
   *                element such that the mapping function is of the form:
   *                \f$\vec{x}(\vec{\xi})=\sum_{i=1}^n N_i(\vec{\xi})\vec{x_i}\f$.
   *                For each vertex indicated in vertex_indices_out, 
   *                this list contains the partial derivatives of the cooresponding
   *                coefficient function \f$N_i\f$ with respect to each 
   *                component of \f$\vec{\xi}\f$ in the same order as the
   *                corresponding nodes in vertex_indices_out. 
   */
  /*@{*/

  /**Get partial derivatives of mapping function coefficients for an 
   * \f$\vec{\xi}\f$ corresponding
   * to the corner of the element.
   *\param corner The element corner, specified as the position of the
   *              corresponding vertex in the canonical ordering of the
   *              element.
   */
  virtual 
  void derivatives_at_corner( unsigned corner, 
                              unsigned nodebits,
                              msq_std::vector<size_t>& vertex_indices_out,
                              msq_std::vector<double>& d_coeff_d_xi_out,
                              MsqError& err ) const = 0;

  /**Get partial derivatives of mapping function coefficients for an 
   * \f$\vec{\xi}\f$ corresponding
   * to the middle of an edge of the element.
   *\param edge The element edge, specified using the canoncial ordering
   *            of edges for the element topology.
   */
  virtual 
  void derivatives_at_mid_edge( unsigned edge, 
                                unsigned nodebits,
                                msq_std::vector<size_t>& vertex_indices_out,
                                msq_std::vector<double>& d_coeff_d_xi_out,
                                MsqError& err ) const = 0;

  /**Get partial derivatives of mapping function coefficients for an 
   * \f$\vec{\xi}\f$ corresponding
   * to the middle of an face of the element.
   *\param face The element face, specified using the canoncial ordering
   *            of faces for the element topology.
   */
  virtual 
  void derivatives_at_mid_face( unsigned face, 
                                unsigned nodebits,
                                msq_std::vector<size_t>& vertex_indices_out,
                                msq_std::vector<double>& d_coeff_d_xi_out,
                                MsqError& err ) const = 0;

  /**Get partial derivatives of mapping function coefficients for an 
   * \f$\vec{\xi}\f$ corresponding
   * to the center of the element.
   */
  virtual 
  void derivatives_at_mid_elem( unsigned nodebits,
                                msq_std::vector<size_t>& vertex_indices_out,
                                msq_std::vector<double>& d_coeff_d_xi_out,
                                MsqError& err ) const = 0;
  /*@}*/
};

} // namespace Mesquite

#endif
