/*!
  \file   LaplacianSmoother.cpp
  \brief  

  The LaplacianSmoother Class is the concrete class
  that performs Laplacian Smoothing

  \author Thomas Leurent
  \date   2002-01-17
*/

#include "LaplacianSmoother.hpp"
#include "MsqMessage.hpp"
#include "LPtoPTemplate.hpp"
#include "EdgeLengthQualityMetric.hpp"


using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::LaplacianSmoother" 
LaplacianSmoother::LaplacianSmoother(MsqError &err) 
{
  this->set_name("LaplacianSmoother");
  
  set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err,1,1);MSQ_CHKERR(err);

  edgeQM = new EdgeLengthQualityMetric;
  edgeQM->set_averaging_method(QualityMetric::RMS,err);
  objFunc = new LPtoPTemplate(edgeQM, 2, err);
  
}  
#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::~LaplacianSmoother" 
LaplacianSmoother::~LaplacianSmoother() 
{
  delete edgeQM;
  delete objFunc;
}    
  
#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::initialize" 
void LaplacianSmoother::initialize(PatchData& /*pd*/, MsqError& /*err*/)
{
 
}

#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::initialize_mesh_iteration" 
void LaplacianSmoother::initialize_mesh_iteration(PatchData &/*pd*/,
                                                  MsqError &/*err*/)
{
  //  cout << "- Executing LaplacianSmoother::iteration_complete()\n";
}

#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::optimize_vertex_position"
/*! \todo Michael:  optimize_vertex_position is probably not implemented
  in an optimal way.  We used to use all of the vertices in
  the patch as 'adjacent' vertices.  Now we call get_adjacent_vertex_indices.
  We could use a VERTICES_ON_VERTEX type of patch or a global patch?
*/
void LaplacianSmoother::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
    //default the laplacian smoother to 3 even for 2-d elements.
    //int dim = get_mesh_set()->space_dim();
  size_t dim = 3;
  
  
  // does the Laplacian smoothing
  MsqFreeVertexIndexIterator free_iter(&pd, err);
  free_iter.reset();
  free_iter.next();
    //m is the free vertex.
  size_t m=free_iter.value();
  std::vector<size_t> vert_indices;
  vert_indices.reserve(25);
    //get vertices adjacent to vertex m
  pd.get_adjacent_vertex_indices(m,vert_indices,err);
    //move vertex m
  centroid_smooth_mesh(pd, vert_indices.size(), vert_indices,
                       m, dim, err); MSQ_CHKERR(err);
    //snap vertex m to domain
  pd.snap_vertex_to_domain(m,err);
  
}
  
#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::terminate_mesh_iteration" 
void LaplacianSmoother::terminate_mesh_iteration(PatchData &/*pd*/,
                                                 MsqError &/*err*/)
{
  //  cout << "- Executing LaplacianSmoother::iteration_complete()\n";
}
  
#undef __FUNC__
#define __FUNC__ "LaplacianSmoother::cleanup" 
void LaplacianSmoother::cleanup()
{
  //  cout << "- Executing LaplacianSmoother::iteration_end()\n";
}
  

