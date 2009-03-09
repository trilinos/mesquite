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


/** \file CachingTargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "CachingTargetCalculator.hpp"
#include "PatchData.hpp"
#include "MsqMatrix.hpp"
#include "ExtraData.hpp"
#include "ElemSampleQM.hpp"

namespace Mesquite {

CachingTargetCalculator::~CachingTargetCalculator() 
{ }

void CachingTargetCalculator::notify_new_patch( PatchData&, CachedTargetData& data )
  { data.clear(); }

void CachingTargetCalculator::notify_patch_destroyed( CachedTargetData& data )
  { data.clear(); }

void CachingTargetCalculator::notify_sub_patch( PatchData& ,
                                                CachedTargetData& data,
                                                PatchData& subpatch, 
                                                const size_t* ,
                                                const size_t* element_map,
                                                MsqError& err )
{
    // If no cached data for this patch, just return
  if (data.has_data())
    return;
  
    // Create a new cached data object on the subpatch
  CachedTargetData& sub_data = get_data( subpatch );
  sub_data.clear();
  
    // populate the element offset list, and count the total
    // number of cached target matrices.
  sub_data.elementOffsets.resize( subpatch.num_elements() );
  size_t count_2D = 0, count_3D = 0;
  for (size_t i = 0; i < subpatch.num_elements(); ++i) {
    EntityTopology type = subpatch.element_by_index(i).get_element_type();
    size_t& count = (TopologyInfo::dimension( type ) == 2) ? count_2D : count_3D;
    sub_data.elementOffsets[i] = count;
    NodeSet samples = subpatch.get_samples( i );
    count += samples.num_nodes();
  }
  sub_data.targets2D.resize( count_2D );
  sub_data.targets3D.resize( count_3D );

  for (size_t i = 0; i < subpatch.num_elements(); ++i) {
    EntityTopology type = subpatch.element_by_index(i).get_element_type();
    size_t off = sub_data.elementOffsets[i];
    size_t old_off = data.elementOffsets[element_map[i]];
    NodeSet samples = subpatch.get_samples( i );
    size_t count = samples.num_nodes();
   
    if (TopologyInfo::dimension( type ) == 3) 
      memcpy( &sub_data.targets3D[off], &data.targets3D[old_off], count*sizeof(MsqMatrix<3,3>) );
    else
      memcpy( &sub_data.targets2D[off], &data.targets2D[old_off], count*sizeof(MsqMatrix<3,2>) );
  }
}

static void populate_data( PatchData& pd,
                           CachedTargetData* data,
                           TargetCalculator* calc,
                           MsqError& err )
{
  size_t i, j;
  
  if (data->elementOffsets.empty()) {
    size_t count_3d = 0, count_2d = 0;
    data->elementOffsets.resize( pd.num_elements() );
    for (i = 0; i < pd.num_elements(); ++i) {
      EntityTopology type = pd.element_by_index(i).get_element_type();
      NodeSet sample_pts = pd.get_samples( i );
      size_t& count = (TopologyInfo::dimension( type ) == 3) ? count_3d : count_2d;
      data->elementOffsets[i] = count;
      count += sample_pts.num_nodes();
    }
    data->targets2D.resize( count_2d );
    data->targets3D.resize( count_3d );
  }
  
  size_t off = 0;
  for (i = 0; i < pd.num_elements(); ++i) {
    EntityTopology type = pd.element_by_index(i).get_element_type();
    NodeSet sample_pts = pd.get_samples( i );
    assert( off == data->elementOffsets[i] );
    if (TopologyInfo::dimension( type ) == 3) {
      for (j = 0; j < TopologyInfo::corners(type); ++j) {
        if (sample_pts.corner_node(j)) {
          calc->get_3D_target( pd, i, Sample(0,j), data->targets3D[off++], err );
          MSQ_ERRRTN(err);
        }
      }
      for (j = 0; j < TopologyInfo::edges(type); ++j) {
        if (sample_pts.mid_edge_node(j)) {
          calc->get_3D_target( pd, i, Sample(1,j), data->targets3D[off++], err );
          MSQ_ERRRTN(err);
        }
      }
      for (j = 0; j < TopologyInfo::faces(type); ++j) {
        if (sample_pts.mid_face_node(j)) {
          calc->get_3D_target( pd, i, Sample(2,j), data->targets3D[off++], err );
          MSQ_ERRRTN(err);
        }
      }
      if (sample_pts.mid_region_node()) {
        calc->get_3D_target( pd, i, Sample(3,0), data->targets3D[off++], err );
        MSQ_ERRRTN(err);
      }
    }
    else {
      for (j = 0; j < TopologyInfo::corners(type); ++j) {
        if (sample_pts.corner_node(j)) {
          calc->get_2D_target( pd, i, Sample(0,j), data->targets2D[off++], err );
          MSQ_ERRRTN(err);
        }
      }
      for (j = 0; j < TopologyInfo::edges(type); ++j) {
        if (sample_pts.mid_edge_node(j)) {
          calc->get_2D_target( pd, i, Sample(1,j), data->targets2D[off++], err );
          MSQ_ERRRTN(err);
        }
      }
      if (sample_pts.mid_face_node(0)) {
        calc->get_2D_target( pd, i, Sample(2,0), data->targets2D[off++], err );
        MSQ_ERRRTN(err);
      }
    }
  }
}

  
bool CachingTargetCalculator::get_3D_target( PatchData& pd, 
                                             size_t element,
                                             Sample sample,
                                             MsqMatrix<3,3>& W_out,
                                             MsqError& err )
{
  CachedTargetData& data = get_data( pd );
  if (data.targets3D.empty()) {
    populate_data( pd, &data, cachedCalculator, err );
    MSQ_ERRZERO(err);
  }
  
    // calculate index of sample in array 
  NodeSet all_samples = pd.get_samples( element );
  unsigned offset = all_samples.num_before( sample );

  W_out = data.targets3D[ data.elementOffsets[element] + offset ];
  return true;
}
  
bool CachingTargetCalculator::get_2D_target( PatchData& pd, 
                                             size_t element,
                                             Sample sample,
                                             MsqMatrix<3,2>& W_out,
                                             MsqError& err )
{
  CachedTargetData& data = get_data( pd );
  if (data.targets2D.empty()) {
    populate_data( pd, &data, cachedCalculator, err );
    MSQ_ERRZERO(err);
    if (data.targets2D.empty()) {
      MSQ_SETERR(err)("Attempt to get 2D target for 3D element type", MsqError::INVALID_STATE);
      return false;
    }
  }
  
    // calculate index of sample in array 
  NodeSet all_samples = pd.get_samples( element );
  unsigned offset = all_samples.num_before( sample );

  W_out = data.targets2D[ data.elementOffsets[element] + offset ];
  return true;
}

} // namespace Mesquite
