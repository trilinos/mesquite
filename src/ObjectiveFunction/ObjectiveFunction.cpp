/*!
  \file   ObjectiveFunction.cpp
  \brief  

  \author Michael Brewer
  \date   2002-08-02
*/

#include "ObjectiveFunction.hpp"
#include "MsqVertex.hpp"
#include "MsqMessage.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "ObjectiveFunction::compute_numerical_gradient"
/*! \fn ObjectiveFunction::compute_numerical_gradient(Mesquite::PatchData &pd, Vector3D *const &grad, MsqError &err, int array_size)
  
  Numerically Calculates the gradient of the ObjectiveFunction for the
  free vertices in the patch.
  \param pd  PatchData on which the gradient is taken.
  \param grad  Array of Vector3D used to store gradient.
  \param array_size Either the length of grad or -1.
 */
bool ObjectiveFunction::compute_numerical_gradient(Mesquite::PatchData &pd,
                                                   Vector3D *const &grad,
                                                   MsqError &err,
                                                   int array_size) {
  int n=pd.num_free_vertices(err); MSQ_CHKERR(err);
  if(n!=array_size && array_size>0)
    PRINT_ERROR("\nArray size not equal n.\n");
  
  MsqVertex* vertices=pd.get_vertex_array(err);
    //index of current free vertex in pd's vertex array
  int m=0;
    //index into gradient array that corresponds to current free vertex
  int grad_pos=0;
  double flocal=0;
  double flocald=0;
  double eps=0;
  MsqFreeVertexIndexIterator ind(&pd, err);
  ind.reset();
  while(ind.next()){
    m=ind.value();
      //This is NOT right, but it is a temporary solution
      //we need an effective way to evaluate an objective function
      //on a local patch... should we create a local patch data object here
      //
      //If pd is not in the feasible region, do not calculate anything.
      //Just return false.
    if(! evaluate(pd,flocal,err))
    {
      return false;
    }
    MSQ_CHKERR(err);
    int j=0;
      //loop over the three coords x,y,z
    for(j=0;j<3;++j){
      eps=get_eps(pd,flocald, j,(&vertices[m]),err);
        //PRINT_INFO("\nin obj num grad j=%i, eps=%20.19f",j,eps);
      if(eps==0){
        err.set_msg("Dividing by zero in Objective Functions numerical grad");
        return false;
      }
      grad[grad_pos][j]=(flocald-flocal)/eps;
    }
    ++grad_pos;
    MSQ_CHKERR(err);
  }//end loop over free vertices
  return true;
  
}

/*!Returns an appropiate value (eps) to use as a delta step for
  MsqVertex vertex in dimension k (i.e. k=0 -> x, k=1 -> y, k=2 -> z).
  The objective function value at the perturbed vertex position is given
  in local_val.
  */
double ObjectiveFunction::get_eps(PatchData &pd, double &local_val,
                                  int k,MsqVertex* vertex, MsqError &err)
{
  double eps = 1.e-07;
//  double rho=.5;
  int imax=20;
  int i=0;
  bool feasible=false;
  double tmp_var=0.0;
  while (i<imax && !feasible)
  {
    i++;
      //perturb kth coord val and check feas if needed
    tmp_var=(*vertex)[k];
    (*vertex)[k]+=eps;
    feasible = evaluate(pd,local_val,err);
      //if step was too big, shorten it         
    if(!feasible)
      eps*=0.5;
      //revert kth coord val
    (*vertex)[k]=tmp_var;
  }//end while looking for feasible eps
 return eps;
}//end funciton get_eps


