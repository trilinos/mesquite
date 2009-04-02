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


/** \file TargetMetric2D.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_METRIC_2D_HPP
#define MSQ_TARGET_METRIC_2D_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {

class MsqError;
template <unsigned R, unsigned C> class MsqMatrix;

/**\brief A metric for comparing a 2x2 matrix A with a 2x2 target matrix W
 *
 * Implement a scalar function \f$\mu(A,W)\f$ where A and W are 2x2 matrices.
 */
class TargetMetric2D {
public:
    // used by code templatized to work with either this class or TargetMetric2D
  enum { MATRIX_DIM = 2 };

  MESQUITE_EXPORT virtual
  ~TargetMetric2D();

    /**\brief Evaluate \f$\mu(A,W)\f$
     *
     *\param A 2x2 active matrix
     *\param W 2x2 target matrix
     *\param result Output: value of function
     *\return false if function cannot be evaluated for given A and W
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& A, 
                 const MsqMatrix<2,2>& W, 
                 double& result, 
                 MsqError& err ) = 0;
  
    /**\brief Gradient of \f$\mu(A,W)\f$ with respect to components of A
     *
     *\param A 2x2 active matrix
     *\param W 2x2 target matrix
     *\param result Output: value of function
     *\param deriv_wrt_A Output: partial deriviatve of \f$\mu\f$ wrt each term of A,
     *                           evaluated at passed A.
     *                           \f[\left[\begin{array}{cc} 
     *                            \frac{\partial\mu}{\partial A_{0,0}} & 
     *                            \frac{\partial\mu}{\partial A_{0,1}} \\ 
     *                            \frac{\partial\mu}{\partial A_{1,0}} & 
     *                            \frac{\partial\mu}{\partial A_{1,1}} \\ 
     *                            \end{array}\right]\f]
     *\return false if function cannot be evaluated for given A and W
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<2,2>& A,
                           const MsqMatrix<2,2>& W,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_A,
                           MsqError& err );

    /**\brief Hessian of \f$\mu(A,W)\f$ with respect to components of A
     *
     *\param A 2x2 active matrix
     *\param W 2x2 target matrix
     *\param result Output: value of function
     *\param deriv_wrt_A Output: partial deriviatve of \f$\mu\f$ wrt each term of A,
     *                           evaluated at passed A.
     *\param second_wrt_A Output: 4x4 matrix of second partial deriviatve of \f$\mu\f$ wrt 
     *                           each term of A, in row-major order.  The symmetric 
     *                           matrix is decomposed into 2x2 blocks and only the upper diagonal
     *                           blocks, in row-major order, are returned.
     *                           \f[\left[\begin{array}{cc|cc}
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{0,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{1,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{1,1}} \\
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{0,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{1,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{1,1}} \\
     *                           \hline & &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{1,1}} \\
     *                           & &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{1,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,1}^2} \\
     *                            \end{array}\right]\f]
     *        
     *\return false if function cannot be evaluated for given A and W
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<2,2>& A,
                           const MsqMatrix<2,2>& W,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_A,
                           MsqMatrix<2,2> second_wrt_A[3],
                           MsqError& err );
};

} // namespace Mesquite

#endif
