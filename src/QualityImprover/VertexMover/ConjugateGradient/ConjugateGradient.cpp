/*!
  \file   ConjugateGradient.cpp
  \brief  

  The Conjugate Gradient class is a concrete vertex mover
  which performs conjugate gradient minimizaiton.

  \author Michael Brewer
  \date   2002-06-19
*/

#include "ConjugateGradient.hpp"
#include<math.h>
#include "ConditionNumberQualityMetric.hpp"
#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
#include "TSTT_C.h"

/* TODO  (Michae) Notes
   What do we do about cull so that it isn't Metric
   dependent... 
   also fix the < or > problem with the Metric.
   Change check_feasible, we don't want to call Condition_number...
   instead we want to check the signed area of the elements.  The
   problem is 2-d elems. in 3-d.
   Figure out the local vs. global patch problem.  We begin with our
   large patch... do we create a local patch to use in grad(...)
   Go through and do clean_up.
*/

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "ConjugateGradient::ConjugateGradient" 
ConjugateGradient::ConjugateGradient(ObjectiveFunction* objective) 
{
  this->set_name("ConjugateGradient");
  //this->set_patch_depth(-1);
  this->set_patch_depth(1);
  objFunc=objective;
  //Michael:: default to global?
  set_step_size_bound(.001);
  set_gradient_bound(.001);
  set_maximum_iteration(6);
}  
  
  
#undef __FUNC__
#define __FUNC__ "ConjugateGradient::initialize" 
void ConjugateGradient::initialize(PatchData &pd, MsqError &err)
{
  PRINT_INFO("\no  Performing Conjugate Gradient optimization.\n");
  arraySize=5;
  fGrad = new Vector3D[ arraySize ];
  pGrad = new Vector3D[ arraySize ];
  fNewGrad = new Vector3D[ arraySize ];
  mCoord = new Vector3D[ arraySize ];
}

#undef __FUNC__
#define __FUNC__ "ConjugateGradient::initialize_mesh_iteration" 
void ConjugateGradient::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
  int n=pd.num_free_vertices();
  if(n>arraySize){
    delete []fGrad;
    delete []pGrad;
    delete []fNewGrad;
    delete []mCoord;
      //Increase number to try to avoid reallocating
    n+=5;
    arraySize=n;
    fGrad = new Vector3D[ n ];
    pGrad = new Vector3D[ n ];
    fNewGrad = new Vector3D[ n ];
    mCoord = new Vector3D[ n ];
  } 
}

#undef __FUNC__
#define __FUNC__ "ConjugateGradient::optimize_vertex_position"
/*!Performs Conjugate gradient minimization on the PatchData, pd.*/
void ConjugateGradient::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err){
  int is=0;
  MeshSet *vertex_mover_mesh=get_mesh_set();
  int num_local_vertices = pd.num_vertices();
  int num_free_vertices=pd.num_free_vertices();
    // gets the array of coordinates for the patch  
  MsqVertex* vertices=pd.get_vertex_array(err);
    //PRINT_INFO("\ny %f\n",vertices[0][1]);
  int ind;
    //This should be done somewhere else
  for(ind=0;ind<num_local_vertices;ind++){
    if(ind<num_free_vertices)
      vertices[ind].set_vertex_flag(MsqVertex::MSQ_FREE_VERTEX);
    else
      vertices[ind].remove_vertex_flag(MsqVertex::MSQ_FREE_VERTEX);
  }
  //Michael hack for geo mesh
#ifdef MESQUITE_GEO         
  if(num_local_vertices-num_free_vertices<6){
    PRINT_INFO("\nSkipping low valence vertex for geo mesh:  vert = %i, free = %i",num_local_vertices,num_free_vertices);
    if(num_local_vertices==1){
      Vector3D hack_vec=vertices[0]->get_vector();
      PRINT_INFO("\nOne vertex patch :  x=%f, y=%f, z=%f",hack_vec[0],hack_vec[1],hack_vec[2]);
      PRINT_INFO("\nNumber of elements = %        i",pd.num_elements());
    }
    
    return;
  }
#endif
    //Michael hack for boundary nodes
  if(num_local_vertices==pd.num_elements()+2)
    return;
  
    //Doing initial feasiblity check if required by metric
  if(objFunc->get_feasible_constraint()){
    if(check_feasible(pd,err))
      err.set_msg("Invalid metric for smoothing a tangled mesh");
  }
    //Michael cull list
    //Michael hard coded for now... fix this
    //But, the culling method needs to be completely re-done anyway
    //For example, calls to evaluate_element????
    //cull_list(pd,free_vertex_list,beta, err);
  int n=pd.num_free_vertices();
  if(n<1){
    PRINT_INFO("\nEmpty free vertex list in ConjugateGradient\n");
  }

  
   
  double f=0;
    
    //Michael, this isn't equivalent to CUBIT because we only want to check
    //the objective function value of the 'bad' elements
  f=objFunc->evaluate(pd,err);
  MSQ_CHKERR(err);
  objFunc->compute_gradient(pd, fGrad , err, n);
  for(ind=0;ind<n;ind++){
    pGrad[ind]=(-fGrad[ind]);
  }
  int j=0;
  int i=0;
  int m=0;
  double alp=MSQ_MAX_CAP;
  bool inner_criterion=inner_criterion_met(*vertex_mover_mesh,err);
  double grad_norm=MSQ_MAX_CAP;
  while ((i<maxIteration && alp>stepBound && grad_norm>normGradientBound)
         && !inner_criterion){
    ++i;
    int k=0;
      //PRINT_INFO("\nbefore y %f\n",vertices[0][1]);
    alp=get_step(pd,f,k,err);
      //PRINT_INFO("\nafter y %f\n",vertices[0][1]);
      //PRINT_INFO("\nAlp initial, alp = %10.8f",alp);
    MSQ_CHKERR(err);
    if(alp==0){
      ++is;
      for (m=0; m<n; m++) {
        pGrad[m]=(-fGrad[m]);
      }
      alp=get_step(pd,f,k,err);
        //PRINT_INFO("\nAlp was zero, alp = %20.18f",alp);
    }
    if(alp!=0){
      j+=k;
      for(m=0;m<n;m++){
        vertices[m] += (alp * pGrad[m]);
      }
      f=objFunc->evaluate(pd,err);
      MSQ_CHKERR(err);
      
      objFunc->compute_gradient(pd, fNewGrad, err, n);
      
      grad_norm=infinity_norm(fNewGrad,n,err);
      
      double s11=0;
      double s12=0;
      double s22=0;
      
      for (m=0; m<n; ++m) {
        s11+=fGrad[m]%fGrad[m];
        s12+=fGrad[m]%fNewGrad[m];
        s22+=fNewGrad[m]%fNewGrad[m];
      }
      
        // Steepest Descent (takes 2-3 times as long as P-R)
        //double bet=0;
      
        // Fletcher-Reeves (takes twice as long as P-R)
        //double bet = s22/s11;

        // Polack-Ribiere        
        double bet = (s22-s12)/s11;
          //PRINT_INFO("\nbet=%10.8f\n",bet);
      for (m=0; m<n; m++) {
        pGrad[m]=(-fNewGrad[m]+(bet*pGrad[m]));
        fGrad[m]=fNewGrad[m];
      }
    }//end if on alp == 0
    inner_criterion=inner_criterion_met(*vertex_mover_mesh,err);MSQ_CHKERR(err);
  }//end while
    //PRINT_INFO("\n2y %f\n",vertices[0][1]);
    //PRINT_INFO("\nConjagate Gradient complete i=%i alp=%4.2e grad_norm=%4.2e   ",i,alp,grad_norm);
}


#undef __FUNC__
#define __FUNC__ "ConjugateGradient::terminate_mesh_iteration" 
void ConjugateGradient::terminate_mesh_iteration(PatchData &pd, MsqError &err)
{
    //  cout << "- Executing ConjugateGradient::iteration_complete()\n";
}


#undef __FUNC__
#define __FUNC__ "ConjugateGradient::cleanup" 
void ConjugateGradient::cleanup()
{
    //  cout << "- Executing ConjugateGradient::iteration_end()\n";
  delete []fGrad;
  delete []pGrad;
  delete []fNewGrad;
  delete []mCoord;
}

//!Computes a distance to move vertices given an initial position and search direction (stored in data member pGrad).
/*!Returns alp, the double which scales the search direction vector
  which when added to the old nodal positions yields the new nodal
  positions.*/
/*!\TODO Michael NOTE:  int &j is only to remain consisitent with CUBIT for an
  initial test.  It can be removed.*/

double ConjugateGradient::get_step(PatchData &pd,double f0,int &j,
                                   MsqError &err){
    //PRINT_INFO("\nEntering get_step  f_x = %10.8f   f_y = %10.8f   f_z = %10.8f\n",pGrad[0][0],pGrad[0][1],pGrad[0][2]);
  const int n = pd.num_free_vertices();
    //iterator for several for statements
  int m=0;
  MsqVertex* vertices=pd.get_vertex_array(err);
    //initial guess for alp
  double alp=1.0;
  int jmax=100;
  double rho=0.5;
    //feasible=1 implies the mesh is not in the feasible region
  int feasible=1;
  int found=0;
    //f and fnew hold the objective function value
  double f=0;
  double fnew=0;
    //Counter to avoid infinitly scaling alp
  j=0;
    //if we must check feasiblility
  for(m=0;m<n;++m)
    mCoord[m]=vertices[m];
  
  if(objFunc->get_feasible_constraint()){
    while (j<jmax && feasible!=0 && alp>MSQ_MIN){
      ++j;
      for (m=0;m<n;++m){
        vertices[m]=(mCoord[m]+ (alp * pGrad[m]));
      }
      
      feasible=check_feasible(pd,err);
      
      MSQ_CHKERR(err);
        //if not feasible, try a smaller alp (take smaller step)
      if(feasible!=0){
        alp*=rho;
      }
    }//end while ...
  }//end if need to check feasible
  else{//if no need for feasiblity, just move using orig. alp
    feasible=0;
    for (m=0;m<n;++m){
      vertices[m]=(mCoord[m])+ (alp * pGrad[m]);
    }
  }
    //if above while ended due to j>=jmax, no valid step was found.
  if(j>=jmax){
    PRINT_INFO("\nFeasible Point Not Found");
    return 0.0;
  }
    //evaluate objective function
  f=objFunc->evaluate(pd,err);
    //PRINT_INFO("\nOrigninal f %f, first new f = %f",f0,f);
    //if new f is larger than original, our step was too large
  if(f>=f0){
    j=0;
    while (j<jmax && found == 0){
      ++j;
        //alp *= (-sqrt(rho));
      alp *= rho;
      for (m=0;m<n;++m){
        vertices[m]=mCoord[m]+ (alp*pGrad[m]);
      }
      f=objFunc->evaluate(pd,err);
      for(m=0;m<n;++m){
        vertices[m][0]=mCoord[m][0];
        vertices[m][1]=mCoord[m][1];
        vertices[m][2]=mCoord[m][2];
      }
      MSQ_CHKERR(err);
        //if our step has now improved the objective function value
      if(f<f0){
        found=1;
          //Michael delete
          //if(alp<=0.0)
          //return 0.0;
      }
    }//   end while j less than jmax
      //PRINT_INFO("\nj = %d found = %d f = %20.18f f0 = %20.18f\n",j,found,f,f0);
      //if above ended because of j>=jmax, take no step
    if(found==0){
        //PRINT_INFO("alp = %10.8f, but returning zero\n",alp);
      alp=0; 
      return alp;
    }
    else{
      return alp;
    }

    j=0;
      //while shrinking the step improves the objFunc value further,
      //scale alp down.  Return alp, when scaling once more would
      //no longer improve the objFunc value.  
    while(j<jmax){
      ++j;
      alp*=rho;
      for(m=0;m<n;++m){
        vertices[m]=mCoord[m] + (alp*pGrad[m]);
        //michael: changed mCoord to mCoord below... check this
      }
      fnew=objFunc->evaluate(pd,err);
      if(fnew<f){
        f=fnew;
      }
      else{
        for(m=0;m<n;++m){
          vertices[m][0]=mCoord[m][0];
          vertices[m][1]=mCoord[m][1];
          vertices[m][2]=mCoord[m][2];
        }
	alp/=rho;
	return alp;
      }
    }
    for(m=0;m<n;++m){
      vertices[m][0]=mCoord[m][0];
      vertices[m][1]=mCoord[m][1];
      vertices[m][2]=mCoord[m][2];
    }
    return alp;
  }
    //else our new f was already smaller than our original
  else{
    j=0;
      //check to see how large of step we can take
    while (j<jmax && found == 0) {
      ++j;
        //scale alp up (rho must be less than 1)
      alp /= rho;
      for(m=0;m<n;++m){
        vertices[m] = mCoord[m] + (alp * pGrad[m]);
      }
      if(objFunc->get_feasible_constraint())
	feasible=check_feasible(pd,err);
      if ( feasible != 0 ){
         alp *= rho;
         for(m=0;m<n;++m){
           vertices[m][0]=mCoord[m][0];
           vertices[m][1]=mCoord[m][1];
           vertices[m][2]=mCoord[m][2];
         }
         return alp;
      }
      fnew=objFunc->evaluate(pd,err);
      if (fnew<f) { 
          f = fnew; 
       }
       else {
         found=1;
         alp *= rho;
       }
    }
    for(m=0;m<n;++m){
      vertices[m][0]=mCoord[m][0];
      vertices[m][1]=mCoord[m][1];
      vertices[m][2]=mCoord[m][2];
    }
    return alp;
  }
}

/*!Quadratic one-dimensional line search.*/
/*double ConjugateGradient::get_step(PatchData &pd,double f0,int &j,
                                   MsqError &err){
  const double CGOLD = 0.3819660;
  const double ZEPS = 1.0e-10;
  int n=pd.num_free_vertices();
  MsqVertex** vertices=pd.get_vertices_array(err);
  double a,b,d,etemp,fb,fu,fv,fw,fx,p,q,r,tol,tol1,tol2,u,v,w,x,xm;
  double e=0.0;
  d=0.0;
  tol=.001;
  int iter, maxiter;
  maxiter=100;
  a=0;
  b=.125;
  int m=0;
  fb=f0-1.0;
  iter=0;
  //find b such that a b 'should' bracket the min
  while (fb<=f0 && iter<maxiter){
    ++iter;
    b*=2.0;
    for(m=0;m<n;++m){
      mCoord[m]=mCoord[m] + (b*pGrad[m]);
      vertices[m]->set_coordinates(mCoord[m]);
    }
    fb=objFunc->evaluate(pd,err);
  }
  iter=0;
  x=w=v=(b/2.0);
  for(m=0;m<n;++m){
    mCoord[m]=mCoord[m] + (x*pGrad[m]);
    vertices[m]->set_coordinates(mCoord[m]);
  }
  fw=fv=fx=objFunc->evaluate(pd,err);
  for(iter=0;iter<maxiter;++iter){
      //PRINT_INFO("a=%f,b=%f,x=%f,iter=%i\n",a,b,x,iter);
    xm=(a+b)*.5;
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if(fabs(x-xm)<= (tol2-0.5*(b-a))){
      return x;
    }
    if(fabs(e)>tol1){
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if(q>0.0)
        p=-p;
      q=fabs(q);
      etemp=e;
      e=d;
      if(fabs(p)>=fabs(0.5*q*etemp)||(p<=q*(a-x))||(p>=q*(b-x))){
        d=CGOLD*(e=(x>=xm?a-x:b-x));
      }
      else{
        d=p/q;
        u=x+d;
        if(u-a<tol2||b-u<tol2)
        {
          if(tol1<0.0)
            d=x-xm;
          else
            d=xm-x;
        }
      }
    }
    
    else{
      d=CGOLD*(e=(x>=xm?a-x:b-x));
    }
    if(tol<0.0)
      u=(fabs(d)>=tol1?x+d:x-d);
    else
      u=(fabs(d)>=tol1?x+d:x+d);
    for(m=0;m<n;++m){
      mCoord[m]=mCoord[m] + (u*pGrad[m]);
      vertices[m]->set_coordinates(mCoord[m]);
    }
    fu=objFunc->evaluate(pd,err);
    if(fu<fx){
      if(u>=x)
        a=x;
      else
        b=x;
      v=w;
      w=x;
      x=u;
      fv=fw;
      fw=fx;
      fx=fu;
    }
    else{
      if(u<x)
        a=u;
      else
        b=u;
      if(fu<=fw||w==x){
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      }
      else if (fu<=fv||v==x||v==w){
        v=u;
        fv=fu;
      }
    }
  }
  for(m=0;m<n;++m){   
    vertices[m]->set_coordinates(mCoord[m]);
  }
  PRINT_WARNING("TOO MANY ITERATIONS IN QUADRATIC LINE SEARCH");
  return x;
}
*/
          
    
/*!cull_list is culls the list of Vertices such that every vertex
  which is attatched to an element with a metric value below beta
  is aded to the list
*/
/*
void ConjugateGradient::cull_list(PatchData &pd,
				double beta, MsqError &err){
  std::list<QualityMetric*> metriclist = objFunc->get_quality_metric_list();
  QualityMetric* qmetric;
  qmetric=metriclist.front();
  MsqMeshEntity* elems = pd.get_elements_array(err);
  int num_elements=pd.num_elements();
  int ind=0;
  for(ind=metriclist.size();ind>0;ind--){
    int jnd=0;
    for(jnd=0;jnd<num_elements;jnd++){
      //TODO (Michael) Fix culling
      //if(qmetric->evaluate_element(elements[jnd],err)<beta){
        list<MsqVertex*> loc_list;
        elems[jnd].get_vertices(loc_list);
        int knd=0;
        for(knd=loc_list.size();knd>0;knd--){
          MsqVertex* tem=loc_list.front();
          loc_list.pop_front();
          if(tem->is_free_vertex()){
              //free_vertex_list.push_back(tem);
          }
        }
	//}
          //free_vertex_list.unique();
    }
  }

}
*/


