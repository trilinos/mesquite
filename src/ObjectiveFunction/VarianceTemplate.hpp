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


/** \file VarianceTemplate.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_VARIANCE_TEMPLATE_HPP
#define MSQ_VARIANCE_TEMPLATE_HPP

#include "Mesquite.hpp"
#include "MsqHessian.hpp"
#include "ObjectiveFunctionTemplate.hpp"

namespace MESQUITE_NS {

/**\brief (standard deviation)^2 template
 *
 * This class implements an objective function that is the 
 * standard deviation of the quality metric evalutations.
 */
class MESQUITE_EXPORT VarianceTemplate : public ObjectiveFunctionTemplate
{
  public:
  
    VarianceTemplate( QualityMetric* qm = 0 ) : ObjectiveFunctionTemplate(qm) 
      { clear(); }
    
      /**\brief copy constructor 
       *
       * Define a copy constructor because the compiler-provided 
       * default one would also copy the temporary arrays, which
       * would be a waste of time.
       */
    VarianceTemplate( const VarianceTemplate& copy )
      : ObjectiveFunctionTemplate( copy ),
        mCount( copy.mCount ),
        mSum( copy.mSum ),
        mSqrSum( copy.mSqrSum ),
        saveCount( copy.saveCount ),
        saveSum( copy.saveSum ),
        saveSqrSum( copy.saveSqrSum )
      {}
    
    virtual ~VarianceTemplate() 
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

    virtual bool evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        msq_std::vector<Vector3D>& grad_out,
                                        msq_std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err ); 

    virtual ObjectiveFunction* clone() const;

    virtual void clear();

  private:
  
    /**\brief Handle EvalType for all eval functions, return OF value 
     *
     * This function implements the common handling of the EvalType
     * argument for all forms of the 'evaluate' method.  
     *
     * NOTE:  This function modifies accumulated values depenending
     *        on the value of EvalType.
     *\param sum        The sum over the current patch
     *\param sqr_sum    The sum of squares over the current patch
     *\param count      The number of qm evaluations for the current patch
     *\param type       The evaluation type passed to 'evaluate'
     *\param global_count The total, accumulated number of QM evaluations
     *\param result_sum The sum term of the standard deviation
     *\param result_sqr The sum of squares term of the standard deviation
     */
    void accumulate( double sum, double sqr_sum, size_t count, 
                     EvalType type,
                     double& result_sum, double& result_sqr, size_t& global_count );
    
    size_t mCount;    /**< The number of accumulated entires */
    double mSum;      /**< The runnnig sum of the qualtiy metric valuse */
    double mSqrSum;   /**< The running sum of the square of QM values */
    size_t saveCount; /**< Saved count from previous patch */
    double saveSum;   /**< Saved sum from previous patch */
    double saveSqrSum;/**< Saved sum from previous patch */
  
    /** Temporary storage for qm sample handles */
    mutable msq_std::vector<size_t> qmHandles;
    /** Temporary storage for qm vertex indices */
    mutable msq_std::vector<size_t> mIndices;
    /** Temporary storage for qm gradient */
    mutable msq_std::vector<Vector3D> mGradient, tmpGradient;
    /** Temporary storage for qm Hessian diagonal data */
    mutable msq_std::vector<SymMatrix3D> mHessDiag, tmpDiag1, tmpDiag2;
};

} // namespace Mesquite

#endif
