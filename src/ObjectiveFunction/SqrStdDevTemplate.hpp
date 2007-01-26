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


/** \file SqrStdDevTemplate.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SQR_STD_DEV_TEMPLATE_HPP
#define MSQ_SQR_STD_DEV_TEMPLATE_HPP

#include "Mesquite.hpp"
#include "StdDevTemplate.hpp"
#include "MsqHessian.hpp"

namespace Mesquite {

/**\brief (standard deviation)^2 template
 *
 * This class implements an objective function that is the 
 * standard deviation of the quality metric evalutations.
 */
class MESQUITE_EXPORT SqrStdDevTemplate : public StdDevTemplate
{
  public:
  
    SqrStdDevTemplate( QualityMetric* qm = 0 ) : StdDevTemplate(qm) {}
    
    SqrStdDevTemplate( const SqrStdDevTemplate& other ) : StdDevTemplate(other) {}
    
    virtual ~SqrStdDevTemplate() 
      {}
    
    virtual bool evaluate( EvalType type, 
                           PatchData& pd,
                           double& value_out,
                           bool free,
                           MsqError& err ); 

    virtual bool evaluate_with_gradient( EvalType type, 
                                         PatchData& pd,
                                         double& value_out,
                                         msq_std::vector<Vector3D>& grad_out,
                                         MsqError& err ); 

    virtual bool evaluate_with_Hessian( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        msq_std::vector<Vector3D>& grad_out,
                                        MsqHessian& Hessian_out,
                                        MsqError& err ); 

    virtual ObjectiveFunction* clone() const;

  private:
  
    /** Temporary storage for qm gradient */
    mutable msq_std::vector<Matrix3D> mHessian;
    mutable MsqHessian tmpHessian1, tmpHessian2;
};

} // namespace Mesquite

#endif
