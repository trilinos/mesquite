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

//#include "TSTT_C.h"

namespace Mesquite
{

  /*! \class VertexMover
    Base class for all Vertex Movers.
   */  
  class VertexMover : public QualityImprover 
  {
  protected:
    VertexMover();
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
    int check_feasible(PatchData &pd, MsqError &err);
  };

  
#undef __FUNC__
#define __FUNC__ "VertexMover::check_feasible"
/*!
  Takes a PatchData object (by reference) and returns whether the
  patch is within the feasible region, 0, or outside the region, 1.
*/
  inline int VertexMover::check_feasible(PatchData &pd, MsqError &err)
  {
    MsqMeshEntity* elems=pd.get_element_array(err);
    int num_elements=pd.num_elements();
    std::vector<Vector3D> sample_points;
    Vector3D jacobian_vectors[3];
    int num_jacobian_vectors;
    int i =0;
    for(i=0;i<num_elements;++i){
      elems[i].get_sample_points(QualityMetric::ELEMENT_VERTICES,sample_points,err);
      std::vector<Vector3D>::iterator iter=sample_points.begin();
      while(iter!=sample_points.end()){
        elems[i].compute_weighted_jacobian(pd, (*iter),
                                           jacobian_vectors,
                                           num_jacobian_vectors, err);
        if(num_jacobian_vectors==2){
            //2-d not yet implemented
        }
        else if(num_jacobian_vectors==3){
          if(jacobian_vectors[0]%(jacobian_vectors[1]*
                                   jacobian_vectors[2])<=0.0){
            return 1;
          }
        }
        ++iter;
      }
    }
    
    return 0;
  }
    
      
        

      
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
        if(grad_norm<fabs(vec[gi][gj])){
          grad_norm=fabs(vec[gi][gj]);
        }
      }
    }
    return grad_norm;
  }

} // namespace
#endif // Mesquite_VertexMover_hpp
