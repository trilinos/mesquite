/*!
  \file   Randomize.hpp
  \brief  

  The Randomize Class implements the Randomize Vertex Mover
  for a patch with one free vertex. 

  \author Michael Brewer      
  \date   2002-10-27
*/

#ifndef Mesquite_Randomize_hpp 
#define Mesquite_Randomize_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "MsqMessage.hpp"

#include <math.h>
namespace Mesquite
{

  /*! \class Randomize
   \brief Randomly perftubs the (un-culled) vertices.
  */ 
  class Randomize : public VertexMover 
  {
  public:
      //!Constructor defaulting mPercent to .05.
    Randomize();
      //!Constructor allowing user to set mPercent
    Randomize(double percent);

  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();
  private:
      //! \param The percentage of the scale factor each vertex will be moved.
    double mPercent;
      /*!Function calculates a scale factor for the patch, then moves
        the incident vertex randomly in each of the three coordinate
        directions (relative to the scale factor multiplied by mPercent).
      */
    void randomize_vertex(int num_incident_vtx,
                          MsqVertex *incident_vtx,
                          MsqVertex &free_vtx,
                          int dimension, MsqError &err);
  };

  
#undef __FUNC__
#define __FUNC__ "randomize_vertex"
    //!Perturbs the free vertex randomly.
  inline void Randomize::randomize_vertex(int num_incident_vtx,
                                          MsqVertex *incident_vtx,
                                          MsqVertex &free_vtx,
                                          int dimension, MsqError &err)
  {
    int i,j;
      //a scale w.r.t. the patch size
    double scale_factor=0.0;
      //a "random" number between -1 and 1
    double rand_double=0.0;
      //a "random" int
    int rand_int=0;
    if (num_incident_vtx==0){
      PRINT_WARNING("WARNING: Number of incident vertex is zero.  Returning.\n");
      return;
    }
    
    for (i=0;i<num_incident_vtx;++i){
      scale_factor+=(incident_vtx[i]-free_vtx).length();
    }
    scale_factor/=( (double) num_incident_vtx );    
    for (j=0;j<3;++j){
      rand_int = rand();
        //number between 0 and 1000
      rand_int = rand_int%1000;
        //number between -1 and 1
      rand_double = (((double) rand_int)/500.0)-1.0;
      free_vtx[j] += scale_factor*rand_double*mPercent;
    }
    
    return;
  }

  
}

#endif
