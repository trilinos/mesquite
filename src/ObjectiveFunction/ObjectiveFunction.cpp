/*!
  \file   ObjectiveFunction.cpp
  \brief  

  \author Michael Brewer
  \author Thomas Leurent
  
  \date   2002-08-02
*/

#include "ObjectiveFunction.hpp"
#include "MsqVertex.hpp"
#include "MsqMessage.hpp"
#include "MsqFreeVertexIndexIterator.hpp"

using namespace Mesquite;
using std::cout;
using std::endl;
using std::cerr;

#undef __FUNC__
#define __FUNC__ "ObjectiveFunction::compute_numerical_gradient"
/*! 
  Numerically Calculates the gradient of the ObjectiveFunction for the
  free vertices in the patch.  Returns 'false' if the patch is outside
  of a required feasible region, returns 'ture' otherwise.
  The behavior of the function depends on the value of the boolean
  useLocalGradient.  If useLocalGradient is set to
  'true', compute_numerical_gradient creates a sub-patch around a free
  vertex, and then perturbs that vertex in one of the coordinate directions.
  Only the ObjectiveFunction value on the local sub-patch is used in the
  computation of the gradient.  Therefore, useLocalGradient should only
  be set to 'true' for ObjectiveFunctions which can use this method.  Unless
  the concrete ObjectiveFunction sets useLocalGradient to 'true' in its
  constructor, the value will be 'false'.  In this case, the objective
  function value for the entire patch is used in the calculation of the
  gradient.  This is computationally expensive, but it is numerically
  correct for all (C_1) functions.
  \param pd  PatchData on which the gradient is taken.
  \param grad  Array of Vector3D of length the number of vertices used to store gradient.
  \param array_size Either the length of grad or 0.
 */
bool ObjectiveFunction::compute_numerical_gradient(Mesquite::PatchData &pd,
                                                   Vector3D *const &grad,
                                                   MsqError &err,
                                                   size_t array_size)
{
  size_t num_vtx=pd.num_vertices();
  if(num_vtx!=array_size && array_size>0)
    PRINT_ERROR("\nArray size not equal to the number of vertices.\n");
  
  MsqVertex* vertices=pd.get_vertex_array(err);
  double flocal=0;
  double flocald=0;
  double eps=0;
  size_t m=0;
  short j;
  
  if(useLocalGradient){
    //********************useLocalGradient***************************
    //if useLocalGradient is turned on, do more efficient computation
    PatchData sub_patch;
    for (m=0; m<num_vtx; ++m) {
     if (vertices[m].is_free_vertex()) {
        pd.get_subpatch(m, sub_patch, err); MSQ_CHKERR(err);
        //If sub_patch is not in the feasible region, do not
        //calculate anything.  Just return false.
        if(! evaluate(sub_patch,flocal,err)) {
          return false;
        }
        MSQ_CHKERR(err);
        //loop over the three coords x,y,z
        for(j=0;j<3;++j){
          eps=get_eps(sub_patch, flocald, j, (&vertices[m]), err);
          //PRINT_INFO("\nin obj num grad j=%i, eps=%20.19f",j,eps);
          if(eps==0){
            err.set_msg("Dividing by zero in Objective Functions numerical grad");
            return false;
          }
          grad[m][j]=(flocald-flocal)/eps;
        }
        MSQ_CHKERR(err);
      }
      else {
        for(j=0;j<3;++j)
          grad[m][j] = 0.0;
      }
    }
  }
  else {
    //********************DO NOT useLocalGradient********************
    //if useLocalGradient is turned off, we do inefficient computation
    for (m=0; m<num_vtx; ++m) {
      if (vertices[m].is_free_vertex()) {
        //If pd is not in the feasible region, do not calculate anything.
        //Just return false.
        if(! evaluate(pd,flocal,err)) {
          return false;
        }
        MSQ_CHKERR(err);
        //loop over the three coords x,y,z
        for(j=0;j<3;++j){
          eps=get_eps(pd, flocald, j, (&vertices[m]), err);
          //PRINT_INFO("\nin obj num grad j=%i, eps=%20.19f",j,eps);
          if(eps==0){
            err.set_msg("Dividing by zero in Objective Functions numerical grad");
            return false;
          }
          grad[m][j]=(flocald-flocal)/eps;
        }
        MSQ_CHKERR(err);
      }
      else {
        for(j=0;j<3;++j)
          grad[m][j] = 0.0;
      }
      //PRINT_INFO("  gradx = %f, grady = %f, gradz = %f\n",grad[m][0],grad[m][1],grad[m][2]);   
      MSQ_CHKERR(err);
    }//end loop over all vertices
  }
  //*****************END of DO NOT useLocalGradient*****************
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
}//end function get_eps


#undef __FUNC__
#define __FUNC__ "ObjectiveFunction::compute_numerical_hessian"
/*! 
  \param pd  PatchData on which the hessian is taken.
  \param  MsqHessian  hessian object. The MsqHessian object needs at least one call to
  MsqHessian::initialize() before being used. 
 */
bool ObjectiveFunction::compute_numerical_hessian(Mesquite::PatchData &/*pd*/,
                                                  MsqHessian &/*hessian*/,
                                                  MsqError &/*err*/)
{
  cout << " THIS FUNCTION EXISTS FOR TEST PURPOSES ONLY.\n"
    "  analytical Objective Function Hessians should always be used. ";
  
//   int num_vtx=pd.num_vertices();
//   Vector3D* grad = new Vector3D[num_vtx];
//   Vector3D* grad_fd = new Vector3D[num_vtx];
//   Vector3D zero(0., 0., 0.);
//   Matrix3D* block;
  
//   this->compute_gradient(pd, grad, err); MSQ_CHKERR(err);
  
//   MsqVertex* vertices=pd.get_vertex_array(err);
//   double flocal=0;
//   double flocald=0;
//   double eps=0;
//   double tmp_coords;
//   int m, v, j;
//   bool grad_success;
  
//   // loop over all vertices in the patch,
//   for (m=0; m<num_vtx; ++m) {
//     //loop over the three coords x,y,z
//     for(j=0;j<3;++j){
//       if (vertices[m].is_free_vertex()) {
//         //        eps=get_eps(pd, flocald, j, (&vertices[m]), err);
//         eps = 10e-6;
//         tmp_coords = vertices[m][j];
//         vertices[m][j] += eps/4;
//         //If pd is not in the feasible region, do not calculate anything.
//         //Just return false.
//         if(eps==0)
//           err.set_msg("Dividing by zero in Objective Function numerical hessian");
//         grad_success = this->compute_gradient(pd, grad_fd, err); MSQ_CHKERR(err);
//         vertices[m][j] = tmp_coords;
//         if( !grad_success ) {
//           cout << "invalid compute_gradient m=" << m << " j=" << j <<endl;//dbg
//           delete[] grad;
//           delete[] grad_fd;
//           return false;
//         }
//         for (v=0; v<=m; ++v) {
//           grad_fd[v] = (grad_fd[v]-grad[v])/eps;
//           cout << "grad_fd["<<v<<"] = "<<grad_fd[v]<<endl;//dbg
//           block = hessian.get_block(v,m);
//           if ( block!=NULL )
//             block->set_column(j, grad_fd[v]);
//         }
//       } else {
//         cout << "vertices["<<m<<"] is fixed.\n";//dbg
//         for (v=m; v<num_vtx; ++v) {
//           hessian.get_block(v,m)->set_column(j, zero);
//         }
//       }
      
//     }
//   }//end loop over all vertices

//   delete[] grad;
//   delete[] grad_fd;
  return true;
}
