
#include "SmartLaplacianSmoother.hpp"
#include "PatchData.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"


namespace Mesquite {

size_t SmartLaplacianSmoother::num_inverted( PatchData& pd, MsqError& err )
{
  MsqMeshEntity::ElementOrientation tmp;
  size_t result = 0;
  for (size_t i = 0; i < pd.num_elements(); ++i) {
    tmp = pd.element_by_index(i).check_element_orientation( pd, err ); 
    MSQ_ERRZERO(err);
    if (tmp != MsqMeshEntity::VALID_ORIENTATION)
      ++result;
  }
  return result;
}


msq_std::string SmartLaplacianSmoother::get_name() const
  { return "SmartLaplacianSmoother"; }

SmartLaplacianSmoother::~SmartLaplacianSmoother() 
{
}    

void SmartLaplacianSmoother::optimize_vertex_positions( PatchData &pd, 
                                                        MsqError &err )
{
  assert(pd.num_free_vertices() == 1);
  const size_t center_vtx_index = 0;
  const size_t init_inverted = num_inverted( pd, err ); MSQ_ERRRTN(err);
  
  adjVtxList.clear();
  pd.get_adjacent_vertex_indices( center_vtx_index, adjVtxList, err );
  MSQ_ERRRTN(err);
  
  if (adjVtxList.empty())
    return;
  
  MsqVertex* verts = pd.get_vertex_array(err);
  const size_t n = adjVtxList.size();
  
  const Vector3D orig_pos = verts[center_vtx_index];
    // static_cast to Vector3D so that we assign only the coorindate data,
    // and not the flags.
  verts[center_vtx_index] = static_cast<Vector3D>(verts[ adjVtxList[0] ]);
  for (size_t i = 1; i < n; ++i)
    verts[center_vtx_index] += verts[ adjVtxList[i] ];
  verts[center_vtx_index] *= 1.0/n;
  
  pd.snap_vertex_to_domain( center_vtx_index, err );  MSQ_ERRRTN(err);
  const size_t new_inverted = num_inverted( pd, err ); MSQ_ERRRTN(err);
  if (new_inverted > init_inverted)
    verts[center_vtx_index] = orig_pos;
}

}
