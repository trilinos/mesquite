/*!
  \file   LaplacianSmoother.hpp
  \brief  

  The LaplacianSmoother Class implements the Laplacian smoothing
  for a patch with one free vertex. 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_LaplacianSmoother_hpp 
#define Mesquite_LaplacianSmoother_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include <vector>
namespace Mesquite
{

  /*! \class LaplacianSmoother
    Moves free center vertex to the average of the neighboring vertices.
   */  
  class LaplacianSmoother : public VertexMover 
  {
  public:
    LaplacianSmoother(MsqError &err);
    ~LaplacianSmoother();
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

  private:
    QualityMetric* edgeQM;

  };

  
#undef __FUNC__
#define __FUNC__ "centroid_smooth_mesh" 
  inline void centroid_smooth_mesh(PatchData &pd, size_t num_adj_vtx,
                                   std::vector<size_t> adj_vtx_ind,
                                   size_t free_ind,
                                   size_t dimension, MsqError &err)
  {
    MsqVertex* verts=pd.get_vertex_array(err);
    std::vector<size_t>::iterator iter;
    
    size_t j;
    double scale_val=1.0;
    if (num_adj_vtx==0) 
      err.set_msg("Number of incident vertices is zero\n");
    else
      scale_val=1.0/((double) num_adj_vtx);
    double avg[3];
      //loop over the two or three dimensions
    for(j=0;j<dimension;++j) {
        //set the iterator to the beginning ob adj_vtx_ind
      iter=adj_vtx_ind.begin();
      avg[j] = 0.;
      while(iter != adj_vtx_ind.end()){
        avg[j]+=verts[*iter][j];
        ++iter;
      }
        //divide the average by the number of adj. verts
      verts[free_ind][j] = avg[j]*scale_val;
    }

    return;
  }

  
}

#endif
