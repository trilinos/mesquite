/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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
 
    (2007) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file AffineMapMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AffineMapMetric.hpp"
#include "MsqMatrix.hpp"
#include "ElementQM.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MappingFunction.hpp"
#include "WeightCalculator.hpp"
#include "TargetCalculator.hpp"
#include "TargetMetric2D.hpp"
#include "TargetMetric3D.hpp"
#include "TargetMetricUtil.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <functional.h>
# include <algorithm.h>
#else
# include <functional>
# include <algorithm>
#endif

namespace Mesquite {

const double TRI_XFORM_VALS[] = { 1.0, -1.0/sqrt(3.0), 0.0, 2.0/sqrt(3.0) };
MsqMatrix<2,2> TRI_XFORM( TRI_XFORM_VALS );

const double TET_XFORM_VALS[] = { 1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0),
                                  0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0),
                                  0.0,  0.0,            sqrt(3.0/2.0) };
MsqMatrix<3,3> TET_XFORM( TET_XFORM_VALS );
 
AffineMapMetric::AffineMapMetric( TargetCalculator* tc,
                                  WeightCalculator* wc,
                                  TargetMetric2D* metric_2d,
                                  TargetMetric3D* metric_3d ) 
  : targetCalc(tc),
    weightCalc(wc),
    metric2D( metric_2d ),
    metric3D( metric_3d )
{ }
 
AffineMapMetric::AffineMapMetric( TargetCalculator* tc,
                                  TargetMetric2D* metric_2d,
                                  TargetMetric3D* metric_3d ) 
  : targetCalc(tc),
    weightCalc(0),
    metric2D( metric_2d ),
    metric3D( metric_3d )
{ }
     

int AffineMapMetric::get_negate_flag( ) const { return 1; }

msq_std::string AffineMapMetric::get_name() const
  { return msq_std::string("AffineMap"); }

void AffineMapMetric::get_evaluations( PatchData& pd,
                                       msq_std::vector<size_t>& handles,
                                       bool free,
                                       MsqError& err )
{
  get_sample_pt_evaluations( pd, handles, free, err );
}

void AffineMapMetric::get_element_evaluations( PatchData& pd,
                                               size_t elem,
                                               msq_std::vector<size_t>& handles,
                                               MsqError& err )
{
  get_elem_sample_points( pd, elem, handles, err );
}

bool AffineMapMetric::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  Sample s = ElemSampleQM::sample( handle );
  size_t e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  unsigned edim = TopologyInfo::dimension( type );
  const size_t* conn = elem.get_vertex_index_array();
  
    // This metric only supports sampling at corners, except for simplices.
    // If element is a simpex, then the Jacobian is constant over a linear 
    // element.  In this case, always evaluate at any vertex.
  unsigned corner = s.number;
  if (s.dimension != 0) {
    if (type == TRIANGLE || type == TETRAHEDRON)
      corner = 0;
    else {
      MSQ_SETERR(err)("Invalid sample point for AffineMapMetric", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  }
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    if (!metric3D) {
      MSQ_SETERR(err)("No 3D metric for affine-map-based metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    Vector3D c[3] = { Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0) };
    unsigned n;
    const unsigned* adj = TopologyInfo::adjacent_vertices( type, s.number, n );
    c[0] = pd.vertex_by_index( conn[adj[0]] ) - pd.vertex_by_index( conn[s.number] );
    c[1] = pd.vertex_by_index( conn[adj[1]] ) - pd.vertex_by_index( conn[s.number] );
    c[2] = pd.vertex_by_index( conn[adj[2]] ) - pd.vertex_by_index( conn[s.number] );
    MsqMatrix<3,3> A( (MsqMatrix<3,1>*)c );
    if (type == TETRAHEDRON)
      A = A * TET_XFORM;

    MsqMatrix<3,3> W;
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    rval = metric3D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
  }
  else {
    if (!metric2D) {
      MSQ_SETERR(err)("No 2D metric for Jacobian-based metric.\n", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  
    Vector3D c[2] = { Vector3D(0,0,0), Vector3D(0,0,0) };
    unsigned n;
    const unsigned* adj = TopologyInfo::adjacent_vertices( type, s.number, n );
    c[0] = pd.vertex_by_index( conn[adj[0]] ) - pd.vertex_by_index( conn[s.number] );
    c[1] = pd.vertex_by_index( conn[adj[1]] ) - pd.vertex_by_index( conn[s.number] );
    MsqMatrix<3,2> App( (MsqMatrix<3,1>*)c );
    
    MsqMatrix<3,2> Wp;
    targetCalc->get_2D_target( pd, e, s, Wp, err ); MSQ_ERRZERO(err);

    MsqMatrix<2,2> A, W;
    MsqMatrix<3,2> RZ;
    surface_to_2d( App, Wp, W, RZ );
    A = transpose(RZ) * App;
    if (type == TRIANGLE)
      A = A * TRI_XFORM;
    
    rval = metric2D->evaluate( A, W, value, err ); MSQ_ERRZERO(err);
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }
  return rval;
}

bool AffineMapMetric::evaluate_with_indices( PatchData& pd,
                                             size_t handle,
                                             double& value,
                                             msq_std::vector<size_t>& indices,
                                             MsqError& err )
{
  Sample   s = ElemSampleQM::sample( handle );
  size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  const size_t* conn = elem.get_vertex_index_array();
  
    // this metric only supports sampling at corners
  if (s.dimension != 0) {
    if (type != TRIANGLE && type != TETRAHEDRON) {
      MSQ_SETERR(err)("Invalid sample point for AffineMapMetric", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    s.dimension = 0;
    s.number = 0;
  }

  unsigned n;
  const unsigned* adj = TopologyInfo::adjacent_vertices( type, s.number, n );
  indices.clear();
  if (conn[s.number] < pd.num_free_vertices())
    indices.push_back(conn[s.number]);
  for (unsigned i = 0; i < n; ++i)
    if (conn[adj[i]] < pd.num_free_vertices())
      indices.push_back(conn[adj[i]]);
  
  return evaluate( pd, handle, value, err );
}

} // namespace Mesquite
