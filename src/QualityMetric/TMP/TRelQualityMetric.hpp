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


/** \file TRelQualityMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TREL_QUALITY_METRIC_HPP
#define MSQ_TREL_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "TMPQualityMetric.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

class TRel2DMetric;
class TRel3DMetric;

/**\brief Compare targets to mapping function Jacobian matrices
 *
 * A quality metric defined using 2D and 3D target metrics,
 * where the active (A) matrix compared to the target by
 * the underlying metrics is the Jacobian matrix of the
 * mapping function at a given sample point.  For surface
 * elements, A is rotated to align the normal with W, such that
 * both matrices can be reduced from 3x2 to 2x2.
 */
class TRelQualityMetric : public TMPQualityMetric
{
public:
  
  /** Used in tests and other templatized code */
  typedef TRel2DMetric Metric2DType;
  typedef TRel3DMetric Metric3DType;

  /**
   *\param tc   The target calculator 
   *\param wc   The weight calculator
   *\param metric_2d Metric to use for surface elements - may be NULL
   *            if mesh contains only volume elements.
   *\param metric_3d Metric to use for volume elements - may be NULL
   *            if mesh contains only surface elements.
   */
  TRelQualityMetric( TargetCalculator* tc,
                     WeightCalculator* wc,
                     TRel2DMetric* metric_2d,
                     TRel3DMetric* metric_3d ) 
    : TMPQualityMetric(tc,wc),
      metric2D( metric_2d ),
      metric3D( metric_3d )
   {}

  /**
   *\param tc   The target calculator 
   *\param metric_2d Metric to use for surface elements - may be NULL
   *            if mesh contains only volume elements.
   *\param metric_3d Metric to use for volume elements - may be NULL
   *            if mesh contains only surface elements.
   */
  TRelQualityMetric( TargetCalculator* tc,
                     TRel2DMetric* metric_2d,
                     TRel3DMetric* metric_3d ) 
    : TMPQualityMetric(tc,0),
      metric2D( metric_2d ),
      metric3D( metric_3d )
   {}
     
  MESQUITE_EXPORT virtual
  std::string get_name() const;
                 
  MESQUITE_EXPORT virtual
  bool evaluate_with_gradient( PatchData& pd,
                               size_t handle,
                               double& value,
                               std::vector<size_t>& indices,
                               std::vector<Vector3D>& gradient,
                               MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian_diagonal( PatchData& pd,
                                       size_t handle,
                                       double& value,
                                       std::vector<size_t>& indices,
                                       std::vector<Vector3D>& gradient,
                                       std::vector<SymMatrix3D>& Hessian_diagonal,
                                       MsqError& err );
                    
  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian( PatchData& pd,
                              size_t handle,
                              double& value,
                              std::vector<size_t>& indices,
                              std::vector<Vector3D>& gradient,
                              std::vector<Matrix3D>& Hessian,
                              MsqError& err );
  
  TRel2DMetric* get_2d_metric() const { return metric2D; }
  TRel3DMetric* get_3d_metric() const { return metric3D; }
  void set_2d_metric( TRel2DMetric* m ) { metric2D = m; }
  void set_3d_metric( TRel3DMetric* m ) { metric3D = m; }

protected:

  MESQUITE_EXPORT virtual
  bool evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 size_t* indices,
                 size_t& num_indices,
                 MsqError& err );

private:

  TRel2DMetric* metric2D;
  TRel3DMetric* metric3D;
};

} // namespace Mesquite

#endif
