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
  inline void centroid_smooth_mesh(PatchData &pd, size_t num_vtx,
                                   MsqVertex *vtx_array,
                                   size_t free_ind,
                                   size_t dimension, MsqError &err)
  {
    size_t i,j;
    double scale_val=1.0;
    if (num_vtx<=1) 
      err.set_msg("WARNING: Number of incident vertex is zero\n");
    else
      scale_val=1.0/((double) num_vtx -1.0);
    double avg[3];
    for (j=0;j<dimension;++j) {
      avg[j] = 0.;
      for (i=0;i<num_vtx;++i){
          //if we are at the free vertex, skip it
        if(i!=free_ind){
          avg[j]+=vtx_array[i][j];
        }
        
      }
      vtx_array[free_ind][j] = avg[j]*scale_val;
    }

    return;
  }

  
}

#endif
