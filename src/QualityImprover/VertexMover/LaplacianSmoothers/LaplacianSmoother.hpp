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
  inline void centroid_smooth_mesh(PatchData &pd, int num_vtx,
                                   MsqVertex *incident_vtx,
                                   MsqVertex &free_vtx,
                                   int dimension, MsqError &err)
  {
    int i,j;
    double avg[3];
    MsqFreeVertexIndexIterator free_iter(&pd, err);
      //figure out which vertex is center
    free_iter.reset();
    free_iter.next();
    int skip_ind=free_iter.value();
    if (num_vtx<=1) 
      err.set_msg("WARNING: Number of incident vertex is zero\n");
      //std::cout << "centroid_smooth_mesh(): original ["<<free_vtx[0]<<","<<free_vtx[1]<<","<<free_vtx[2]<<"] = " << std::endl;
    for (j=0;j<dimension;++j) {
      free_iter.reset();
      free_iter.next();
      skip_ind=free_iter.value();
      avg[j] = 0.;
      for (i=0;i<num_vtx;++i){
          //cout<<"INSIDE v:  i="<<i<<"   "<<incident_vtx[i];
          //if we are at the free vertex, skip it
        if(i==skip_ind){
          free_iter.next();
          skip_ind=free_iter.value();
        }
          //otherwise:
        else{
          avg[j]+=incident_vtx[i][j];
        }
      }
      free_vtx[j] = avg[j]/((double) num_vtx - 1.0);
        //std::cout << "centroid_smooth_mesh(): final --  avg["<<j<<"] = " << free_vtx[j] << std::endl;
    }

    return;
  }

  
}

#endif
