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

  edgeQM = EdgeLengthQualityMetric::create_new();
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
void LaplacianSmoother::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
    //std::cout << "- Executing LaplacianSmoother::optimize_vertex_position()\n";
    //if(pd.num_free_vertices(err)!=1){
    //PRINT_INFO("\nNum free vertices equals %i\n",pd.num_free_vertices(err));
    //}
  
  int num_local_vertices = pd.num_vertices();
    //default the laplacian smoother to 3 even for 2-d elements.
    //int dim = get_mesh_set()->space_dim();
  size_t dim = 3;
  
  
  
  // gets the array of coordinates for the patch and print it 
  MsqVertex *patch_coords = pd.get_vertex_array(err); MSQ_CHKERR(err);
    //for (size_t i=0; i<num_local_vertices; i++) 
      //cout << "vertex " << i << " : " << patch_coords[i];
  // does the dumb Laplacian smoothing
  MsqFreeVertexIndexIterator free_iter(&pd, err);
  free_iter.reset();
  free_iter.next();
    //find the free vertex.
  size_t m=free_iter.value();
  centroid_smooth_mesh(pd, num_local_vertices, &patch_coords[0],
                       m, dim, err); MSQ_CHKERR(err);
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
  

