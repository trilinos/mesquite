/*!
  \file   VertexMover.hpp
  \brief  

  The VertexMover Class is the base class for all the smoothing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_VertexMover_hpp 
#define Mesquite_VertexMover_hpp


#include "Mesquite.hpp"
#include "QualityImprover.hpp"
#include "PatchData.hpp"

#include "TSTT_C.h"

namespace Mesquite
{

  /*   */ 
  class VertexMover : public QualityImprover 
  {
  public:
    // virtual destructor ensures use of polymorphism during destruction
    virtual ~VertexMover() { };
    
    virtual void loop_over_mesh(MeshSet &ms, MsqError &err);

  protected:

    virtual void initialize(PatchData &pd, MsqError &err) = 0;
    virtual void cleanup() = 0;
    virtual void optimize_vertex_positions(PatchData &pd, 
                                           MsqError &err) = 0; // modifies the PatchData object

    virtual void initialize_mesh_iteration(PatchData &pd, 
                                         MsqError &err) = 0;
    virtual void terminate_mesh_iteration(PatchData &, 
                                         MsqError &err) = 0;
    //!Computes the L_inf norm of an array of Vector3D of length len
    double infinity_norm(Vector3D * const vec, int len, MsqError &err);

      //!CHECK FEASIBLE IS NOT YET IMPLEMENTED.
    int check_feasible(PatchData &pd, MsqError &err)
      {
        return 0;
      }
    


  };
#undef __FUNC__
#define __FUNC__ "VertexMover::infinity_norm"
/*!
  Takes an array of Vecort3D, vec, of length len, and
  and returns the infinity norm of vec.
*/
  inline double VertexMover:: infinity_norm(Vector3D * const vec, int len,
                                          MsqError &err)
  {
    double grad_norm=0;
    for(int gi=0;gi<len;++gi){
      for (int gj=0;gj<3;++gj){
        if(grad_norm<abs(vec[gi][gj])){
          grad_norm=abs(vec[gi][gj]);
        }
      }
    }
    return grad_norm;
  }

} // namespace
#endif // Mesquite_VertexMover_hpp
