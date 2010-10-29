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


/** \file TRel2DMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TREL_2D_METRIC_HPP
#define MSQ_TREL_2D_METRIC_HPP

#include "Mesquite.hpp"
#include <string>

namespace MESQUITE_NS {

class MsqError;
template <unsigned R, unsigned C> class MsqMatrix;

/**\brief A metric for evaluating the relative quality of a matrix as compared
 *        to a desired target metrix
 *
 * Implement a scalar function \f$\mu(T)\f$ where T is a 2x2 matrix.
 */
class TRel2DMetric {
public:
    // used by code templatized to work with either this class or TRel2DMetric
  enum { MATRIX_DIM = 2 };

  MESQUITE_EXPORT virtual
  ~TRel2DMetric();

  MESQUITE_EXPORT virtual
  std::string get_name() const = 0;

    /**\brief Evaluate \f$\mu(T)\f$
     *
     *\param T 2x2 relative measure matrix (typically A W^-1)
     *\param result Output: value of function
     *\return false if function cannot be evaluated for given T
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& T, 
                 double& result, 
                 MsqError& err ) = 0;
  
    /**\brief Gradient of \f$\mu(T)\f$ with respect to components of T
     *
     *\param T 2x2 relative measure matrix (typically A W^-1)
     *\param result Output: value of function
     *\param deriv_wrt_T Output: partial deriviatve of \f$\mu\f$ wrt each term of T,
     *                           evaluated at passed T.
     *                           \f[\left[\begin{array}{cc} 
     *                            \frac{\partial\mu}{\partial T_{0,0}} & 
     *                            \frac{\partial\mu}{\partial T_{0,1}} \\ 
     *                            \frac{\partial\mu}{\partial T_{1,0}} & 
     *                            \frac{\partial\mu}{\partial T_{1,1}} \\ 
     *                            \end{array}\right]\f]
     *\return false if function cannot be evaluated for given T
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<2,2>& T,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_T,
                           MsqError& err );

    /**\brief Hessian of \f$\mu(T)\f$ with respect to components of T
     *
     *\param T 3x3 relative measure matrix (typically A W^-1)
     *\param result Output: value of function
     *\param deriv_wrt_T Output: partial deriviatve of \f$\mu\f$ wrt each term of T,
     *                           evaluated at passed T.
     *\param second_wrt_T Output: 9x9 matrix of second partial deriviatve of \f$\mu\f$ wrt 
     *                           each term of T, in row-major order.  The symmetric 
     *                           matrix is decomposed into 3x3 blocks and only the upper diagonal
     *                           blocks, in row-major order, are returned.
     *                           \f[\left[\begin{array}{cc|cc}
     *                           \frac{\partial^{2}\mu}{\partial T_{0,0}^2} &
     *                           \frac{\partial^{2}\mu}{\partial T_{0,0}\partial A_{0,1}} &
     *                           \frac{\partial^{2}\mu}{\partial T_{0,0}\partial A_{1,0}} &
     *                           \frac{\partial^{2}\mu}{\partial T_{0,0}\partial A_{1,1}} \\
     *                           \frac{\partial^{2}\mu}{\partial T_{0,0}\partial A_{0,1}} &
     *                           \frac{\partial^{2}\mu}{\partial T_{0,1}^2} &
     *                           \frac{\partial^{2}\mu}{\partial T_{0,1}\partial A_{1,0}} &
     *                           \frac{\partial^{2}\mu}{\partial T_{0,1}\partial A_{1,1}} \\
     *                           \hline & &
     *                           \frac{\partial^{2}\mu}{\partial T_{1,0}^2} &
     *                           \frac{\partial^{2}\mu}{\partial T_{1,0}\partial A_{1,1}} \\
     *                           & &
     *                           \frac{\partial^{2}\mu}{\partial T_{1,0}\partial A_{1,1}} &
     *                           \frac{\partial^{2}\mu}{\partial T_{1,1}^2} \\
     *                            \end{array}\right]\f]
     *        
     *\return false if function cannot be evaluated for given T
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<2,2>& T,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_T,
                           MsqMatrix<2,2> second_wrt_T[3],
                           MsqError& err );
                           
protected:
  static inline bool invalid_determinant( double d )
    { return d < 1e-12; }
};

} // namespace Mesquite

#endif
