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


/** \file TargetSurfaceOrientation.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_SURFACE_ORIENTATION_HPP
#define MSQ_TARGET_SURFACE_ORIENTATION_HPP

#include "Mesquite.hpp"
#include "ElemSampleQM.hpp"

namespace Mesquite {

class TargetCalculator;
class WeightCalculator;
class SamplePoints;
template <unsigned R, unsigned C> class MsqMatrix;

class TargetSurfaceOrientation : public ElemSampleQM
{
public:

  TargetSurfaceOrientation( 
                  const SamplePoints* pts, 
                  TargetCalculator* tc,
                  WeightCalculator* wc = 0 ) 
    : samplePts(pts), 
      targetCalc(tc),
      weightCalc(wc)
   {}
     
  MESQUITE_EXPORT virtual
  msq_std::string get_name() const;
  
  MESQUITE_EXPORT virtual 
  int get_negate_flag() const;
  
  MESQUITE_EXPORT virtual
  void get_evaluations( PatchData& pd, 
                        msq_std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );
  
  MESQUITE_EXPORT virtual 
  void get_element_evaluations( PatchData& pd, size_t elem_index,
                                msq_std::vector<size_t>& handles,
                                MsqError& err );
  
  MESQUITE_EXPORT virtual
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 msq_std::vector<size_t>& indices,
                 MsqError& err );
                 
  MESQUITE_EXPORT 
  const SamplePoints* get_sample_points() const 
    { return samplePts; }
  
private:
  const SamplePoints* samplePts;
  TargetCalculator* targetCalc;
  WeightCalculator* weightCalc;
  
  msq_std::vector<size_t> mIndices;
  msq_std::vector<double> mDerivs;
};

} // namespace Mesquite

#endif
