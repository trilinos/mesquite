/*!
  \file   ObjectiveFunction.cpp
  \brief  

  \author Michael Brewer
  \date   2002-08-02
*/

#include "ObjectiveFunction.hpp"
#include "MsqVertex.hpp"
#include "MsqMessage.hpp"
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
void ObjectiveFunction::compute_numerical_gradient(Mesquite::PatchData &pd,
                                                   Vector3D *const &grad,
                                                   MsqError &err,
                                                   int array_size) {
  int n=pd.num_free_vertices(err); MSQ_CHKERR(err);
  if(n!=array_size && array_size>0)
    PRINT_ERROR("\nArray size not equal n.\n");
  
  MsqVertex* vertices=pd.get_vertex_array(err);
  int m=0;
  double flocal=0;
  double flocald=0;
  double eps=0;
  double coord_hold;
//  double xtmp,ytmp,ztmp;
  for(m=0; m<n; m++){
      //This is NOT right, but it is a temporary solution
      //we need an effective way to evaluate an objective function
      //on a local patch... should we create a local patch data object here
    flocal=evaluate(pd,err);
    MSQ_CHKERR(err);
    int j=0;
      //loop over the three coords x,y,z
    for(j=0;j<3;++j){
      eps=get_eps(pd,j,(&vertices[m]),err);
        //PRINT_INFO("\nin obj num grad j=%i, eps=%20.19f",j,eps);
      coord_hold=vertices[m][j];
      vertices[m][j]+=eps;
      flocald=evaluate(pd,err);
        //PRINT_INFO("\nin obj num grad j=%i, flocal=%18.16f, flocald=%18.16f",j,flocal,flocald);
      if(eps==0){
        err.set_msg("Dividing by zero in Objective Functions numerical grad");
        return;
      }
      grad[m][j]=(flocald-flocal)/eps;
      vertices[m][j]=coord_hold;
    }
    MSQ_CHKERR(err);
  }//end loop over free vertices
}

/*!Returns an appropiate value (eps) to use as a delta step for
  MsqVertex vertex in dimension k (i.e. k=0 -> x, k=1 -> y, k=2 -> z).
  */
double ObjectiveFunction::get_eps(PatchData &pd,int k,MsqVertex* vertex,
                                  MsqError &err)
{
  double eps = 1.e-07;
//  double rho=.5;
  int imax=20;
  int i=0;
  int feasible=1;
  double tmp_var=0.0;
  while (i<imax&&feasible==1)
  {
    i++;
      //perturb kth coord val and check feas if needed
    tmp_var=(*vertex)[k];
    (*vertex)[k]+=eps;
    if(get_feasible_constraint()){
      feasible=check_feasible(pd,err);
        //if step was too big, shorten it
      if(feasible==1)
        eps*=0.5;
    }
      //revert kth coord val
    (*vertex)[k]=tmp_var;
  }//end while looking for feasible eps
  return eps;
}//end funciton get_eps


//!Check feasible takes a patch and makes sure that each element is valid
int ObjectiveFunction::check_feasible(PatchData &pd,MsqError &err){
  //TODO (Michael) THIS FUNCTION NEEDS TO BE IMPLEMENTED
  return 0;    
}
