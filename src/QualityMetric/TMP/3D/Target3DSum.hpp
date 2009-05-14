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


/** \file Target3DSum.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET2D_SUM_HPP
#define MSQ_TARGET2D_SUM_HPP

#include "Mesquite.hpp"
#include "TargetMetric3D.hpp"

namespace MESQUITE_NS {

/** \f$ \mu\prime = \mu_1 + \mu_2 \f$ */
class Target3DSum : public TargetMetric3D
{
  TargetMetric3D *mu1, *mu2;

public:

  Target3DSum( TargetMetric3D* metric1, TargetMetric3D* metric2 ) 
    : mu1(metric1), mu2(metric2) {}
  
  MESQUITE_EXPORT virtual
  msq_std::string get_name() const;

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& A, 
                 const MsqMatrix<3,3>& W, 
                 double& result, 
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<3,3>& A,
                           const MsqMatrix<3,3>& W,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_A,
                           MsqError& err );
  
  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<3,3>& A,
                           const MsqMatrix<3,3>& W,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_A,
                           MsqMatrix<3,3> second_wrt_A[6],
                           MsqError& err );
};


} // namespace MESQUITE_NS

#endif
