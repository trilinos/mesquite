/*!  \file NonSmoothSteepestDescent.cpp \brief
  
  Implements the NonSmoothSteepestDescent class member functions.
  
  \author Lori Freitag
  \date 2002-07-20 */

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "NonSmoothSteepestDescent.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::NonSmoothSteepestDescent" 
NonSmoothSteepestDescent::NonSmoothSteepestDescent(ObjectiveFunction* of)
  : objFunc(of)
{
  this->set_name("NonSmoothSteepestDescent");

  function = (double *)malloc(sizeof(double)*150);
  test_function = (double *)malloc(sizeof(double)*150);
  original_function = (double *)malloc(sizeof(double)*150);
  gradient = (double **)malloc(sizeof(double *)*150);
  gs = (double *)malloc(sizeof(double)*150);
  G = (double **)malloc(sizeof(double *)*150);
  PDG = (double **)malloc(sizeof(double *)*3);
  prev_active_values=(double *)malloc(sizeof(double)*150);

  for (int i=0;i<150;i++) {
    gradient[i] = (double *)malloc(sizeof(double)*3);
    G[i] = (double *)malloc(sizeof(double)*150);
  }
  for (int i=0;i<3;i++) PDG[i] = (double *)malloc(sizeof(double)*3);

  active = (ActiveSet *)malloc(sizeof(ActiveSet));
  test_active = (ActiveSet *)malloc(sizeof(ActiveSet));
  original_active = (ActiveSet *)malloc(sizeof(ActiveSet));
  std::cout << "- Executed NonSmoothSteepestDescent::NonSmoothSteepestDescent()\n";
}  
  
  
#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::initialize" 
void NonSmoothSteepestDescent::initialize(PatchData &pd, MsqError &err)
{
  this->set_patch_depth(1);
  
    // local parameter initialization
  max_iterations = 100;
  conv_eps = 1e-10;
  active_epsilon = .00003;
  min_acceptable_improvement = 1e-6;
  min_step_size = 1e-6;
  std::cout << "- Executed NonSmoothSteepestDescent::initialize()\n";
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::initialize_mesh_iteration" 
void NonSmoothSteepestDescent::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::optimize_vertex_positions" 
void NonSmoothSteepestDescent::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  std::cout << "- Executing NonSmoothSteepestDescent::optimize_node_positions()\n";
  /* perform the min max smoothing algorithm */
  MSQ_DEBUG_PRINT(2,"\nInitializing the patch iteration\n");

  numVertices = pd.num_vertices();
  MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Number of Vertices: %d\n",numVertices);});
  numElements = pd.num_elements();
  MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Number of Elements: %d\n",numElements);});
  dimension = get_mesh_set()->space_dim();
  MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Spatial Dimension: %d\n",dimension);});
  coords = pd.get_vertex_array(err); MSQ_CHKERR(err);
  MSQ_DEBUG_ACTION(3,{
    for (int i99=0;i99<numVertices;i99++) {
      fprintf(stdout,"coords: %g %g\n",coords[i99][0],coords[i99][1]);
    }
  });

  connectivity = pd.get_element_array(err); MSQ_CHKERR(err);
  MSQ_DEBUG_ACTION(3,{
    std::vector<size_t> indices;
    for (int i99=0;i99<numElements;i99++) {
      connectivity[i99].get_vertex_indices(indices);
      fprintf(stdout,"connectivity: %d %d %d\n",indices[0],
	      indices[1],indices[2]);
    }
  });

  /* check for an invalid mesh; if it's invalid return and ask the user 
     to use untangle */
  if (this->validity_check(err)!=1) {
      err.set_msg("Invalid Mesh: Use untangle to create a valid triangulation");
      return;
  }

  /* assumes one function value per element */
  num_function_values = numElements;

  /* initialize the optimization data up to num_function_values */
  this->init_opt(err);
  this->init_max_step_length(err); MSQ_CHKERR(err);
  MSQ_DEBUG_PRINT(3,"Done initializing optimization\n");

  /* compute the initial function values */
  this->compute_function(&pd, original_function, err); MSQ_CHKERR(err);
 
  // find the initial active set
  this->find_active_set(original_function, active, err);  MSQ_CHKERR(err);

  this->minmax_opt(pd,err); MSQ_CHKERR(err);
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::terminate_mesh_iteration" 
void NonSmoothSteepestDescent::terminate_mesh_iteration(PatchData &pd, MsqError &err)
{
}
  
#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::cleanup" 
void NonSmoothSteepestDescent::cleanup()
{
  std::cout << "- Executing NonSmoothSteepestDescent::cleanup()\n";
  for (int i=0;i<150;i++) {
    free(gradient[i]);
    free(G[i]);
  }
  for (int i=0;i<3;i++) free(PDG[i]);

  free(function);
  free(test_function);
  free(original_function);
  free(gradient);
  free(gs);
  free(G);
  free(PDG);
  free(prev_active_values);
  free(active);
  free(test_active);
  free(original_active);
  std::cout << "- Done with NonSmoothSteepestDescent::cleanup()\n";
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::compute_function"
void NonSmoothSteepestDescent::compute_function(PatchData *patch_data, double *func, MsqError &err)
{
  // ASSUMES ONE VALUE PER ELEMENT; ALSO NEED 1.0/FUNCTION WHICH IS ONLY
  // TRUE OF CONDITION NUMBER

  //  MSQ_DEBUG_PRINT(2,"Computing Function\n");

  for (int i=0;i<numElements;i++) func[i]=0.0;
  QualityMetric* currentQM=objFunc->get_quality_metric();
  if(currentQM==NULL){
    currentQM = objFunc->get_quality_metric_list().front();
  }
  
  for (int i=0;i<numElements;i++) {
    func[i] = currentQM->evaluate_element(*patch_data,
                                          &(patch_data->element_by_index(i)),
                                          err); MSQ_CHKERR(err);
    func[i] = 1.0/func[i];
    //    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Function value[%d]=%g\n",i,func[i]);});
  }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::compute_gradient"
double** NonSmoothSteepestDescent::compute_gradient(PatchData *patch_data, MsqError &err)
{
  MSQ_DEBUG_PRINT(2,"Computing Gradient\n");

  double delta = 10e-6;

  for (int i=0;i<numElements;i++) {
    for (int j=0;j<3;j++) gradient[i][j] = 0.0;
  }
  QualityMetric* currentQM=objFunc->get_quality_metric();
  if(currentQM==NULL)
    currentQM = objFunc->get_quality_metric_list().front();

  double *func, *fdelta;
  func = (double *)malloc(sizeof(double)*150);
  fdelta = (double *)malloc(sizeof(double)*150);

  this->compute_function(patch_data, func, err); MSQ_CHKERR(err);

  /* gradient in the x, y, z direction */
  for (int j=0;j<3;j++) {

    // perturb the coordinates of the free vertex in the j direction by delta
    coords[0][j] += delta;
      //patch_data->set_coords_array_element(coords[0],0,err);

    //compute the function at the perturbed point location
    this->compute_function(patch_data, fdelta, err); MSQ_CHKERR(err);

    //compute the numerical gradient
    for (int i=0;i<num_function_values;i++) {
       gradient[i][j] = (fdelta[i] - func[i])/delta;
       // MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Gradient value[%d][%d]=%g\n",i,j,gradient[i][j]);});
    }

    // put the coordinates back where they belong
    coords[0][j] -= delta;
      // patch_data->set_coords_array_element(coords[0],0,err);
  }
  
  free(func);
  free(fdelta);
  return(gradient);
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::find_active_set"
void NonSmoothSteepestDescent::find_active_set(double *function, ActiveSet *active_set,
                                               MsqError &err)
{ 
    int         i, ind;
    double      function_val;
    double      active_value0;
    double      temp;

    MSQ_DEBUG_PRINT(2,"\nFinding the active set\n");

    // initialize the active set indices to zero
    for (i=0;i<num_function_values;i++) active_set->active_ind[i] = 0; 

    /* the first function value is our initial active value */
    active_set->num_active = 1;
    active_set->num_equal = 0;
    active_set->active_ind[0] = 0;
    active_set->true_active_value = function[0];
    //    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  function value[0]: %g\n",function[0]);});

    /* first sort out the active set... 
       all vals within active_epsilon of largest val */

    for (i=1;i<num_function_values;i++) {
	function_val = function[i];
        active_set->true_active_value = MSQ_MAX(function_val,active_set->true_active_value);
	active_value0 = function[active_set->active_ind[0]];
	temp = fabs(function_val - active_value0);
	//        MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  function value[%d]: %g\n",i,function[i]);});
	if ( function_val > active_value0 ) {
	    if ( temp > active_epsilon) {
		active_set->num_active = 1;
		active_set->num_equal = 0;
		active_set->active_ind[0] = i;
	    } else if ( temp < active_epsilon) {
		active_set->num_active += 1;
		ind = active_set->num_active - 1;
		active_set->active_ind[ind] = i;
		if (fabs(function_val - active_value0) < MSQ_MACHINE_EPS) {
		    active_set->num_equal += 1;
		}
	    }
	} else {
	    if (temp < active_epsilon) {
		active_set->num_active += 1;
		ind = active_set->num_active - 1;
		active_set->active_ind[ind] = i;
		if (fabs(function_val - active_value0) < MSQ_MACHINE_EPS) {
		    active_set->num_equal += 1;
		}
	    }
	}
    }

    MSQ_DEBUG_ACTION(3,{
      /* Print the active set */
      this->print_active_set(active_set,function,err); MSQ_CHKERR(err);
    });

}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::improvement_check"
int NonSmoothSteepestDescent::improvement_check(MsqError &err)
{
  int valid = 1;
  
  /* check to see that the mesh didn't get worse */
  if (original_value < active->true_active_value) {
       printf("The local mesh got worse; initial value %f; final value %f\n",
	       original_value,  active->true_active_value );
       valid = 0;
   }

  return(valid);

}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::validity_check"
int NonSmoothSteepestDescent::validity_check(MsqError &err)
        
{
  // ONLY FOR SIMPLICIAL MESHES - THERE SHOULD BE A VALIDITY CHECKER ASSOCIATED
  // WITH MSQ ELEMENTS
  
  /* check that the simplicial mesh is still valid, based on right handedness. 
       Returns a 1 or a 0 */
  int valid = 1;
  double dEps = 1.e-13;

  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

  //  MSQ_DEBUG_PRINT(2,"\nChecking Mesh Validity\n");

  if (dimension == 2)
  {
    for (int i=0;i<numElements;i++)
    {
      double dummy = 0;
      coords[connectivity[i].get_vertex_index(0)].get_coordinates(x1, y1, dummy);
      coords[connectivity[i].get_vertex_index(1)].get_coordinates(x2, y2, dummy);
      coords[connectivity[i].get_vertex_index(2)].get_coordinates(x3, y3, dummy);
      //      MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Element %d\n",i);});
      //      MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  x1 y1 %g %g \n",x1,y1);});
      //      MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  x2 y2 %g %g \n",x2,y2);});
      //      MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  x3 y3 %g %g \n",x3,y3);});
      
      double a = x2*y3 - x3*y2;
      double b = y2 - y3;
      double c = x3 - x2;
      
      if (.5*(a+b*x1+c*y1) < .01*MSQ_MACHINE_EPS)
        valid=0;
    }
  }

  if (dimension == 3)
  {
    for (int i=0;i<numElements;i++)
    {
      x1=coords[0][0];     
      y1=coords[0][1];
      z1=coords[0][2];
      
      coords[connectivity[i].get_vertex_index(0)].get_coordinates(x2, y2, z2);
      coords[connectivity[i].get_vertex_index(1)].get_coordinates(x3, y3, z3);
      coords[connectivity[i].get_vertex_index(2)].get_coordinates(x4, y4, z4);
      
      double dDX2 = x2 - x1;
      double dDX3 = x3 - x1;
      double dDX4 = x4 - x1;
      
      double dDY2 = y2 - y1;
      double dDY3 = y3 - y1;
      double dDY4 = y4 - y1;

      double dDZ2 = z2 - z1;
      double dDZ3 = z3 - z1;
      double dDZ4 = z4 - z1;
      
        /* dDet is proportional to the cell volume */
      double dDet = dDX2*dDY3*dDZ4 + dDX3*dDY4*dDZ2 + dDX4*dDY2*dDZ3
        - dDZ2*dDY3*dDX4 - dDZ3*dDY4*dDX2 - dDZ4*dDY2*dDX3 ;

        /* Compute a length scale based on edge lengths. */
      double dScale = ( sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) +
                             (z1-z2)*(z1-z2)) +
                        sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) +
                             (z1-z3)*(z1-z3)) +
                        sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) +
                             (z1-z4)*(z1-z4)) +
                        sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) +
                             (z2-z3)*(z2-z3)) +
                        sqrt((x2-x4)*(x2-x4) + (y2-y4)*(y2-y4) +
                             (z2-z4)*(z2-z4)) +
                        sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) +
                             (z3-z4)*(z3-z4)) ) / 6.;
      
        /* Use the length scale to get a better idea if the tet is flat or
           just really small. */
      if (fabs(dScale) < MSQ_MACHINE_EPS)
      {
        return(valid = 0);
      }
      else
      {
        dDet /= (dScale*dScale*dScale);
      }
      
      if (dDet > dEps)
      {
        valid = 1;
      }
      else if (dDet < -dEps)
      {
        valid = -1;
      }
      else
      {
        valid = 0;
      }
    }  // end for i=1,numElements
  }  // end dimension==3
  
  //  MSQ_DEBUG_ACTION(2,{fprintf(stdout,"Mesh Validity is: %d \n",valid);});
  
  return(valid);
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::check_equilibrium"
void NonSmoothSteepestDescent::check_equilibrium(int *equil, int *status, MsqError &err)
{
//    int  ierr;
    int  i,j;
    int  ind1, ind2;  
    double min;
    double **dir;
    double mid_vec[3], mid_cos, test_cos;

    *equil = MSQ_FALSE;
    ind1 = ind2 = -1;

    int num_active = active->num_active;

    if (num_active==num_function_values)
    {
         *equil = 1; *status = MSQ_EQUILIBRIUM;
         MSQ_DEBUG_PRINT(3,"All the function values are in the active set\n"); 
    }

    /* set up the active gradient directions */
    this->get_active_directions(gradient,&dir,err);

    /* normalize the active directions */
    for (j=0;j<num_active;j++) MSQ_NORMALIZE(dir[j],dimension);

    if (dimension == 2) {
      /* form the grammian */
      this->form_grammian(dir,err);  

    /* find the minimum element in the upper triangular portion of G*/
    min = 1;
    for (i=0;i<num_active;i++) {
      for (j=i+1;j<num_active;j++) {
        if ( fabs(-1 - G[i][j]) < 1E-08 ) {
           *equil = 1; *status = MSQ_EQUILIBRIUM;
           MSQ_DEBUG_PRINT(3,"The gradients are antiparallel, eq. pt\n"); 
         }
         if (G[i][j]  < min) {
           ind1 = i; ind2 = j;
           min = G[i][j];
        }
      }
    }

    if ((ind1 != -1) && (ind2 != -1)) {
      /* find the diagonal of the parallelepiped */
      for (j=0;j<dimension;j++) {
       mid_vec[j]=.5*(dir[ind1][j]+dir[ind2][j]);
      }
      MSQ_NORMALIZE(mid_vec,dimension);
      MSQ_DOT(mid_cos,dir[ind1],mid_vec,dimension);

      /* test the other vectors to be sure they lie in the cone */
      for (i=0;i<num_active;i++) {
         if ((i != ind1) && (i != ind2)) {
            MSQ_DOT(test_cos,dir[i],mid_vec,dimension);
            if ((test_cos < mid_cos)  &&  fabs(test_cos-mid_cos) > MSQ_MACHINE_EPS) {
              MSQ_DEBUG_PRINT(3,"An equilibrium point \n");
              *equil = 1; *status = MSQ_EQUILIBRIUM;
            }
         }
       }
     }
    }
    if (dimension == 3) {
       *equil = this->convex_hull_test(dir,num_active,err);
       if (*equil == 1) *status = MSQ_EQUILIBRIUM;
    }
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::convex_hull_test"
int NonSmoothSteepestDescent::convex_hull_test(double **vec, int num_vec, MsqError &err)
{
//    int ierr;
    int equil;
    int status, dir_done;
    double pt1[3], pt2[3], pt3[3];
    double normal[3];

    /* tries to determine equilibrium for the 3D case */
    equil = 0;
    status = MSQ_CHECK_Z_COORD_DIRECTION;
    dir_done = -1;

    if (num_vec <= 2) status = MSQ_NO_EQUIL;

    while ((status != MSQ_EQUIL) && (status != MSQ_NO_EQUIL) && 
           (status != MSQ_HULL_TEST_ERROR)) {
       if (status == MSQ_CHECK_Z_COORD_DIRECTION) {
          this->find_plane_points(MSQ_ZDIR, MSQ_YDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status, err);
          dir_done = 2;
       }
       if (status == MSQ_CHECK_Y_COORD_DIRECTION) {
          this->find_plane_points(MSQ_YDIR, MSQ_XDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status, err);
          dir_done = 1;
       }
       if (status == MSQ_CHECK_X_COORD_DIRECTION) {
          this->find_plane_points(MSQ_XDIR, MSQ_ZDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status, err);
          dir_done = 0;
       }
       if (status == MSQ_TWO_PT_PLANE) {
          pt3[0]=0.; pt3[1]=0.; pt3[2]=0.;
       }
       if ((status == MSQ_TWO_PT_PLANE) || (status == MSQ_THREE_PT_PLANE)){
           this->find_plane_normal(pt1,pt2,pt3,normal,err); 
           equil = this->check_vector_dots(vec,num_vec,normal,err); 
           if (equil == 1) {
             switch(dir_done){
             case 2:
               equil = 0; status = MSQ_CHECK_Y_COORD_DIRECTION;
               break;
             case 1:
               equil = 0; status = MSQ_CHECK_X_COORD_DIRECTION;
               break;
             case 0:
               equil = 1; status = MSQ_EQUIL;
             }
           } else if (equil == 0) {
               status = MSQ_NO_EQUIL;
           } else {
               err.set_msg("equil flag not set to 0 or 1");
           }
       }
    }
    switch (status){
    case MSQ_NO_EQUIL:
      MSQ_DEBUG_PRINT(3,"Not an equilibrium point\n");
      equil = 0; 
      break;
    case MSQ_EQUIL:
      MSQ_DEBUG_PRINT(3,"An equilibrium point\n");
      equil = 1;
      break;
    default:
      MSQ_DEBUG_ACTION(3,{
        fprintf(stdout,"Failed to determine equil or not; status = %d\n",status);
      });
    }
    return (equil);
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::check_vector_dots"
int NonSmoothSteepestDescent::check_vector_dots(double **vec, int num_vec, 
                                             double *normal, MsqError &err)
{
    int equil;
    int i, ind;
    double test_dot, dot;

    equil = MSQ_FALSE;
    MSQ_DOT(test_dot,vec[0],normal,3);
    ind = 1;
    while ((fabs(test_dot) < MSQ_MACHINE_EPS) && (ind<num_vec)) {
      MSQ_DOT(test_dot,vec[ind],normal,3);
      ind++;
    }
      
    for (i=ind;i<num_vec;i++) {
       MSQ_DOT(dot,vec[i],normal,3);
       if ( ((dot>0 && test_dot<0) || (dot<0 && test_dot>0)) &&
            (fabs(dot)>MSQ_MACHINE_EPS)) {
          return(equil = MSQ_TRUE);

       }
    }
    return(equil);
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::find_plane_normal"
void NonSmoothSteepestDescent::find_plane_normal(double pt1[3], double pt2[3], double pt3[3], 
                       double *cross, MsqError &err)
{
  int i;
  double vecA[3], vecB[3];

  for (i=0;i<3;i++) {
    vecA[i] = pt2[i] - pt1[i];
    vecB[i] = pt3[i] - pt1[i];
  }
  cross[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
  cross[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
  cross[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];
  MSQ_NORMALIZE(cross, 3);
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::find_plane_points"
void NonSmoothSteepestDescent::find_plane_points(int dir1, int dir2, double **vec, 
                          int num_vec, double *pt1,
			  double *pt2, double*pt3, int *status, MsqError &err)
{
    int i;
    int ind[50], num_min, num_max;
    int rotate=MSQ_CW;
    int num_rotated=0;
    double pt_1, pt_2;
    double min, inv_slope;
    double min_inv_slope=0.;
    double max; 
    double max_inv_slope=0;
    double inv_origin_slope=0;

    *status = MSQ_CHECK_BOTTOM_UP;
    /* find the minimum points in dir1 starting at -1 */
    num_min = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; min=1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1]<min) {
	min = vec[i][dir1]; ind[0] = i; num_min = 1;
      } else if (fabs(vec[i][dir1] - min) < MSQ_MACHINE_EPS) {
	ind[num_min++] = i;
      }
    }
    if (min >= 0) *status = MSQ_NO_EQUIL;
 
    if (*status != MSQ_NO_EQUIL) {
      switch(num_min) {
      case 1: /* rotate to find the next point */
	MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	pt_1 = pt1[dir1]; pt_2 = pt1[dir2];
	if (pt1[dir2] <= 0){rotate=MSQ_CCW; max_inv_slope=MSQ_BIG_NEG_NMBR;}
	if (pt1[dir2] > 0){rotate=MSQ_CW; min_inv_slope=MSQ_BIG_POS_NMBR;}
	switch(rotate) {
	case MSQ_CCW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope>max_inv_slope) &&  
		  (fabs(inv_slope - max_inv_slope) > MSQ_MACHINE_EPS)) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case MSQ_CW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope<min_inv_slope) && 
		  (fabs(inv_slope - max_inv_slope) > MSQ_MACHINE_EPS)){
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  MSQ_DEBUG_PRINT(3,"No points in the rotation ... odd\n");
	    *status = MSQ_HULL_TEST_ERROR;
	  break;
	case 1:
	  MSQ_DEBUG_PRINT(3,"Found a line in the convex hull\n");
	  MSQ_COPY_VECTOR(pt2,vec[ind[1]],3); *status = MSQ_TWO_PT_PLANE;
	  break;
	default:
	  MSQ_DEBUG_PRINT(3,"Found 2 or more points in the rotation\n");
	    if (fabs(pt_1) > MSQ_MACHINE_EPS) inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case MSQ_CCW:
	    if (inv_origin_slope >= max_inv_slope) *status=MSQ_NO_EQUIL;
	    else *status=MSQ_CHECK_TOP_DOWN;
	    break;
	  case MSQ_CW:
	    if (inv_origin_slope <= min_inv_slope) *status=MSQ_NO_EQUIL;
	    else *status=MSQ_CHECK_TOP_DOWN;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
	MSQ_DEBUG_PRINT(3,"Found two minimum points to define the plane\n");
                MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	MSQ_COPY_VECTOR(pt2,vec[ind[1]],3);
	*status = MSQ_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Found 3 or more points in min plane %f\n",min);})
	  if (vec[ind[0]][dir1] >= 0) *status = MSQ_NO_EQUIL;
	  else *status = MSQ_CHECK_TOP_DOWN;
    }
    }

    /***************************/
    /*  failed to find any information, checking top/down this coord*/
    /***************************/

    if (*status == MSQ_CHECK_TOP_DOWN) {
    /* find the maximum points in dir1 starting at 1 */
    num_max = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; max=-1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1] > max) {
	max = vec[i][dir1]; ind[0] = i; num_max = 1;
      } else if (fabs(vec[i][dir1] - max) < MSQ_MACHINE_EPS) {
	ind[num_max++] = i;
      }
    }
    if (max <= 0) *status = MSQ_NO_EQUIL;
 
    if (*status != MSQ_NO_EQUIL) {
      switch(num_max) {
      case 1: /* rotate to find the next point */
	MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	pt_1 = pt1[dir1];  pt_2 = pt1[dir2];
	if (pt1[dir2] < 0){rotate=MSQ_CW; min_inv_slope=1E300;}
	if (pt1[dir2] >= 0){rotate=MSQ_CCW; max_inv_slope=-1E300;}
	switch(rotate) {
	case MSQ_CCW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope>max_inv_slope) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case MSQ_CW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope<min_inv_slope) {
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < MSQ_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  MSQ_DEBUG_PRINT(3,"No points in the rotation ... odd\n");
	  *status = MSQ_HULL_TEST_ERROR;
	  break;
	case 1:
	  MSQ_DEBUG_PRINT(3,"Found a line in the convex hull\n");
          MSQ_COPY_VECTOR(pt2,vec[ind[1]],3);
	  *status = MSQ_TWO_PT_PLANE;
	  break;
	default:
	  MSQ_DEBUG_PRINT(3,"Found 2 or more points in the rotation\n");
	    /* check to see if rotation got past origin */
	  inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case MSQ_CCW:
	    if (inv_origin_slope >= max_inv_slope) *status=MSQ_NO_EQUIL;
	    else if (dir1 == 2) *status=MSQ_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) *status=MSQ_CHECK_X_COORD_DIRECTION;
	    else *status=MSQ_EQUIL;
	    break;
	  case MSQ_CW:
	    if (inv_origin_slope <= min_inv_slope) *status=MSQ_NO_EQUIL;
	    else if (dir1 == 2) *status=MSQ_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) *status=MSQ_CHECK_X_COORD_DIRECTION;
	    else *status=MSQ_EQUIL;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
	MSQ_COPY_VECTOR(pt1,vec[ind[0]],3);
	MSQ_COPY_VECTOR(pt2,vec[ind[1]],3);
	*status = MSQ_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Found 3 in max plane %f\n",max);});
	if (vec[ind[0]][dir1] <= 0) *status = MSQ_NO_EQUIL;
	else if (dir1==2) *status=MSQ_CHECK_Y_COORD_DIRECTION;
	else if (dir1==1) *status=MSQ_CHECK_X_COORD_DIRECTION;
	else *status = MSQ_EQUIL;
      }
    }
  }

}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::form_grammian" 
void NonSmoothSteepestDescent::form_grammian(double **vec, MsqError &err)
{
   int i, j;
   int num_active = active->num_active;

   if (num_active > 150) {
      err.set_msg("Exceeded maximum allowed active values");
   }
   /* form the grammian with the dot products of the gradients */
   for (i=0; i<num_active; i++) {
      for (j=i; j<num_active; j++) {
         G[i][j] = 0.;
	 MSQ_DOT(G[i][j],vec[i],vec[j],dimension);
	 G[j][i] = G[i][j];	 
      }
   }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::form_PD_grammian" 
void NonSmoothSteepestDescent::form_PD_grammian(MsqError &err)
{
    int  i,j,k;
    int  g_ind_1;
    int  singular;

    int num_active = active->num_active;
        
    /* this assumes the grammian has been formed */
    for (i=0;i<num_active;i++) {
      for (j=0;j<num_active;j++) {
        if (G[i][j]==-1) {
          err.set_msg("Grammian not computed properly");
          return;
        }
      }
    }

    /* use the first gradient in the active set */
    g_ind_1 = 0;
    PDG[0][0] = G[0][0];
    PDG_ind[0] = active->active_ind[0];

    /* test the rest and add them as appropriate */
    k = 1; i = 1;
    while( (k<dimension) && (i < num_active) ) {
        PDG[0][k] = PDG[k][0] = G[0][i];
        PDG[k][k] = G[i][i];
        if ( k == 2) { /* add the dot product of g1 and g2 */
           PDG[1][k] = PDG[k][1] = G[g_ind_1][i];
        }
        this->singular_test(k+1,PDG,&singular,err);
        if (!singular) {
           PDG_ind[k] = active->active_ind[i];
           if (k==1) g_ind_1 = i;
           k++;
        }
        i++;
    }
    num_LI = k;
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::singular_test" 
void NonSmoothSteepestDescent::singular_test(int n, double **A, int *singular, MsqError &err) 
{
//    int test;
//    double determinant;
    double cond;

    if ((n>3) || (n<1)) {
      err.set_msg("Singular test works only for n=1 to n=3");
    }

    (*singular)=MSQ_TRUE;
    switch(n) {
    case 1:
        if (A[0][0] > 0) (*singular) = MSQ_FALSE;
        break;
    case 2:
        if (fabs(A[0][0]*A[1][1] - A[0][1]*A[1][0]) > MSQ_MACHINE_EPS)           
            (*singular) = MSQ_FALSE;
        break;
    case 3:
       /* calculate the condition number */
        this->condition3x3(A, &cond, err); 
        if (cond < 1E14) (*singular)=MSQ_FALSE;
        break;
    }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::condition3x3" 
void NonSmoothSteepestDescent::condition3x3(double **A, double *cond, MsqError &err) 
{
//   int ierr;
   double a11, a12, a13;
   double a21, a22, a23;
   double a31, a32, a33;
//   double s1, s2, s4, s3, t0;
   double s1, s2, s3;
   double denom;
//   double one = 1.0;
   double temp;
   int zero_denom = MSQ_TRUE;

   a11 = A[0][0]; a12=A[0][1]; a13=A[0][2];
   a21 = A[1][0]; a22=A[1][1]; a23=A[1][2];
   a31 = A[2][0]; a32=A[2][1]; a33=A[2][2];

   denom = -a11*a22*a33+a11*a23*a32+a21*a12*a33-a21*a13*a32-
            a31*a12*a23+a31*a13*a22;

   if ( (fabs(a11) > MSQ_MACHINE_EPS) && 
        (fabs(denom/a11) > MSQ_MACHINE_EPS)) {
         zero_denom = MSQ_FALSE;
   }
   if ( (fabs(a22) > MSQ_MACHINE_EPS) && 
        (fabs(denom/a22) > MSQ_MACHINE_EPS)) {
         zero_denom = MSQ_FALSE;
   }       
   if ( (fabs(a33) > MSQ_MACHINE_EPS) && 
        (fabs(denom/a33) > MSQ_MACHINE_EPS)) {
         zero_denom = MSQ_FALSE;
   }

   if (zero_denom) {
     (*cond) = 1E300;
   } else {
     s1 = sqrt(a11*a11 + a12*a12 + a13*a13 + 
               a21*a21 + a22*a22 + a23*a23 + 
               a31*a31 + a32*a32 + a33*a33);

     temp = (-a22*a33+a23*a32)/denom;
     s3 = temp*temp;
     temp =(a12*a33-a13*a32)/denom;
     s3 += temp*temp;
     temp = (a12*a23-a13*a22)/denom;
     s3 += temp*temp;
     temp = (a21*a33-a23*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a33-a13*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a23-a13*a21)/denom;
     s3 += temp*temp;
     temp = (a21*a32-a22*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a32+a12*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a22+a12*a21)/denom;
     s3 += temp*temp;

     s2 = sqrt(s3);
     (*cond) = s1*s2;
   }
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::get_active_directions" 
void NonSmoothSteepestDescent::get_active_directions(double **gradient, 
                                         double ***dir, MsqError &err)
{
    int i;
    int num_active = active->num_active;

    (*dir) =(double **)malloc(sizeof(double *)*num_active);
    for (i=0;i<num_active;i++) {
        (*dir)[i] =(double *)malloc(sizeof(double)*dimension);
        MSQ_COPY_VECTOR((*dir)[i],gradient[active->active_ind[i]],dimension);
    }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::search_direction" 
void NonSmoothSteepestDescent::search_direction(PatchData &pd, MsqError &err)
{
   int        i;
   int        viable;
   int        singular;
   double     a, b, c, denom;
   double     **dir;
   double     R0, R1;
   double     **P, *x;
   double     search_mag;

   int num_active = active->num_active;

   MSQ_DEBUG_PRINT(2,"\nIn Search Direction\n");
   this->print_active_set(active, function, err);
   
   if (num_active==0) 
       err.set_msg("No active values in search");

    switch(num_active) {
    case 1: 
        MSQ_COPY_VECTOR(search,gradient[active->active_ind[0]],dimension);
        steepest = active->active_ind[0];
        break;
    case 2:
        /* if there are two active points, move in the direction of the
	   intersection of the planes.  This is the steepest descent
           direction found by analytically solving the QP */
        
        /* set up the active gradient directions */
        this->get_active_directions(gradient,&dir,err);  MSQ_CHKERR(err);

        /* form the grammian */
        this->form_grammian(dir,err);  MSQ_CHKERR(err);
        this->form_PD_grammian(err); MSQ_CHKERR(err);

        denom = (G[0][0] + G[1][1] - 2*G[0][1]);
        viable = 1;
        if (fabs(denom) > MSQ_MACHINE_EPS) {
	  /* gradients are LI, move along their intersection */
           b = (G[0][0] - G[0][1])/denom;  
           a = 1 - b;
           if ((b < 0) || (b > 1)) viable=0;  /* 0 < b < 1 */
           if (viable) {
             for (i=0;i<dimension;i++) {
               search[i] = a*dir[0][i] + b*dir[1][i];
             }
           } else {
             /* the gradients are dependent, move along one face */
             MSQ_COPY_VECTOR(search,dir[0],dimension);
           }
        } else {
	   /* the gradients are dependent, move along one face */
           MSQ_COPY_VECTOR(search,dir[0],dimension);
        }
        steepest = active->active_ind[0];

        for (i=0;i<num_active;i++) free(dir[i]);
	free(dir);

        break;
    default:
        /* as in case 2: solve the QP problem to find the steepest
           descent direction.  This can be done analytically - as
           is done in Gill, Murray and Wright 
             for 3 active points in 3 directions - test PD of G
             otherwise we know it's SP SD so search edges and faces */

        /* get the active gradient directions */
        this->get_active_directions(gradient,&dir,err);  MSQ_CHKERR(err);

        /* form the entries of the grammian matrix */
        this->form_grammian(dir,err);  MSQ_CHKERR(err);
        this->form_PD_grammian(err); MSQ_CHKERR(err);

        switch(dimension) {
        case 2:
  	    this->search_edges_faces(dir,err); MSQ_CHKERR(err);
            break;
        case 3:
	  if (num_active == 3) {
              this->singular_test(num_active,G,&singular,err); MSQ_CHKERR(err);
              if (!singular) {
	        /* form the entries of P=Z^T G Z where Z = [-1...-1; I ] */
                this->form_reduced_matrix(&P,err); MSQ_CHKERR(err);
                /* form  the RHS and solve the system for the coeffs */
                R0 = G[0][0] - G[1][0];  R1 = G[0][0] - G[2][0];
                this->solve2x2(P[0][0],P[0][1],P[1][0],P[1][1],R0,R1,&x,err);
                if (x!=NULL) {
                	a = 1 - x[0] - x[1];  b = x[0];  c = x[1];
                	for (i=0;i<dimension;i++) {
                    	  search[i] = a*dir[0][i] + b*dir[1][i] + 
                       	                      c*dir[2][i];
                	}
                        steepest = active->active_ind[0];
                	for (i=0;i<num_active-1;i++)  free(P[i]);  
                	free(P);  free(x);
                } else { 
                  	this->search_edges_faces(dir, err);
                        MSQ_CHKERR(err);
                }
	      } else {
                 this->search_edges_faces(dir, err); MSQ_CHKERR(err);
	      }
            } else {
              this->search_edges_faces(dir, err); MSQ_CHKERR(err);
            }
            break;
        }
        for (i=0;i<num_active;i++) free(dir[i]);
	free(dir);
    }

    /* if the search direction is essentially zero, equilibrium pt */
    MSQ_DOT(search_mag,search,search,dimension);
    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Search Magnitude %g \n",search_mag);});

    if (fabs(search_mag)<1E-13) opt_status = MSQ_ZERO_SEARCH;
    else MSQ_NORMALIZE(search,dimension);
    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Search Direction %g %g  Steepest %d\n",search[0],search[1],steepest);});
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::search_edges_faces" 
void NonSmoothSteepestDescent::search_edges_faces(double **dir, MsqError &err)
{
//    int ierr;
    int i,j,k;
    int viable;
    double a,b,denom;
    double g_bar[3];
    double temp_search[3];
    double projection, min_projection;

    int num_active = active->num_active;

    if ( (dimension != 2) && (dimension != 3)) {
       err.set_msg("Dimension must be 2 or 3");
    }

    /* initialize the search direction to 0,0 */
    for (i=0;i<dimension;i++) temp_search[i] = 0;

    /* Check for viable faces */
    min_projection = 1E300;
    for (i=0; i<num_active; i++) {
        /* FACE I */
        viable = 1;

        /* test the viability */
        for (j=0;j<num_active;j++) {       /* lagrange multipliers>0 */
             if (G[j][i] < 0) viable = 0;
        }
       
        /* find the minimum of viable directions */
        if ((viable) && (G[i][i] < min_projection)) {
            min_projection = G[i][i];
            MSQ_COPY_VECTOR(temp_search,dir[i],dimension);
            steepest = active->active_ind[i];
        }
    
       /* INTERSECTION IJ */
       for (j=i+1; j<num_active; j++) {
          viable = 1;

          /* find the coefficients of the intersection 
             and test the viability */
          denom = 2*G[i][j] - G[i][i] - G[j][j];
          a = b = 0;
          if (fabs(denom) > MSQ_MACHINE_EPS) {
             b = (G[i][j] - G[i][i])/denom;
             a = 1 - b;
             if ((b < 0) || (b > 1)) viable=0;  /* 0 < b < 1 */
	     for (k=0;k<num_active;k++) {       /* lagrange multipliers>0 */
                 if ((a*G[k][i] + b*G[k][j]) <= 0) viable=0;
             }
          } else {
             viable = 0;                        /* Linearly dependent */
          }

          /* find the minimum of viable directions */
          if (viable) {
             for (k=0;k<dimension;k++) {
                g_bar[k] = a * dir[i][k] + b * dir[j][k];
             }
             MSQ_DOT(projection,g_bar,g_bar,dimension);
             if (projection < min_projection) {
	        min_projection = projection;
                MSQ_COPY_VECTOR(temp_search,g_bar,dimension);
                steepest = active->active_ind[i];
             }
          }
       }
    }
    if (opt_status != MSQ_EQUILIBRIUM) {
        MSQ_COPY_VECTOR(search,temp_search,dimension);
    }
}         

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::solve2x2" 
void NonSmoothSteepestDescent::solve2x2(double a11, double a12, double a21, double a22, 
		   double b1, double b2, double **x, MsqError &err)
{
    double factor;

    /* if the system is not singular, solve it */
    if (fabs(a11*a22 - a21*a12) > MSQ_MACHINE_EPS) {
	(*x)=(double *)malloc(sizeof(double)*2);
	if (fabs(a11) > MSQ_MACHINE_EPS) {
	    factor = (a21/a11);
	    (*x)[1] = (b2 - factor*b1)/(a22 - factor*a12);
	    (*x)[0] = (b1 - a12*(*x)[1])/a11;
	} else if (fabs(a21) > MSQ_MACHINE_EPS) {
	    factor = (a11/a21);
	    (*x)[1] = (b1 - factor*b2)/(a12 - factor*a22);
	    (*x)[0] = (b2 - a22*(*x)[1])/a21;
	}
    } else {
	(*x) = NULL;
    }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::form_reduced_matrix" 
void NonSmoothSteepestDescent::form_reduced_matrix(double ***P, MsqError &err)
{
    int i,j;
    int num_active = active->num_active;

    (*P)=(double **)malloc(sizeof(double *)*(num_active-1));
    for (i=0; i<num_active-1; i++) 
        (*P)[i]=(double *)malloc(sizeof(double)*(num_active-1));

    for (i=0;i<num_active-1;i++) {
        (*P)[i][i] = G[0][0] - 2*G[0][i+1] + G[i+1][i+1];
        for (j=i+1;j<num_active-1;j++) {
            (*P)[i][j] = G[0][0] - G[0][j+1] - G[i+1][0] + G[i+1][j+1];
            (*P)[j][i] = (*P)[i][j];
        }
    }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::minmax_opt"
void NonSmoothSteepestDescent::minmax_opt(PatchData &pd, MsqError &err)
{
//      int valid;

      MSQ_DEBUG_PRINT(2,"In minmax_opt\n");

      MSQ_COPY_VECTOR(function,original_function,num_function_values);
      original_value = active->true_active_value;

      iter_count = 0;
      opt_iter_count = 0;

      MSQ_DEBUG_PRINT(3,"Done copying original function to function\n");

      this->find_active_set(function, active, err); MSQ_CHKERR(err);
      prev_active_values[0] = active->true_active_value;

     /* check for equilibrium point */
     /* compute the gradient */
     gradient = this->compute_gradient(&pd, err); MSQ_CHKERR(err);
     
     if (active->num_active >= 2) {
	MSQ_DEBUG_PRINT(3,"Testing for an equilibrium point \n");
	this->check_equilibrium(&equilibrium_pt, &opt_status, err); MSQ_CHKERR(err);

	MSQ_DEBUG_ACTION(2,{
	    if (equilibrium_pt) 
		fprintf(stdout,"Optimization Exiting: An equilibrium point \n");
        });
     }

    /* terminate if we have found an equilibrium point or if the step is
       too small to be worthwhile continuing */
    while ((opt_status != MSQ_EQUILIBRIUM) && 
	   (opt_status != MSQ_STEP_TOO_SMALL) &&
	   (opt_status != MSQ_IMP_TOO_SMALL) &&
	   (opt_status != MSQ_FLAT_NO_IMP) &&
           (opt_status != MSQ_ZERO_SEARCH) &&
	   (opt_status != MSQ_MAX_ITER_EXCEEDED)) {

	/* increase the iteration count by one */
        /* smooth_param->iter_count += 1; */
        iter_count += 1;
        opt_iter_count += 1;
        if (iter_count > MSQ_MAX_OPT_ITER) opt_status = MSQ_MAX_ITER_EXCEEDED;

	MSQ_DEBUG_PRINT(3,"\n");
	MSQ_DEBUG_ACTION(3,{ 
            fprintf(stdout,"ITERATION %d \n",iter_count);
        });
	    
	/* compute the gradient */
	gradient = this->compute_gradient(&pd, err); MSQ_CHKERR(err);
        
	MSQ_DEBUG_PRINT(3,"Computing the search direction \n");
	this->search_direction(pd, err); MSQ_CHKERR(err);

	/* if there are viable directions to search */
	if ((opt_status != MSQ_ZERO_SEARCH) &&
            (opt_status != MSQ_MAX_ITER_EXCEEDED)) {

	    MSQ_DEBUG_PRINT(3,"Computing the projections of the gradients \n");
	    this->get_gradient_projections(err); MSQ_CHKERR(err);

	    MSQ_DEBUG_PRINT(3,"Computing the initial step size \n");
	    this->compute_alpha(err); MSQ_CHKERR(err);

	    MSQ_DEBUG_PRINT(3,"Testing whether to accept this step \n");
	    this->step_acceptance(pd, err); MSQ_CHKERR(err);
            MSQ_DEBUG_ACTION(3,
              {printf("The new free vertex position is %f %f %f\n",
              coords[0][0],coords[0][1],coords[0][2]);});

	    MSQ_DEBUG_ACTION(3,{
     		/* Print the active set */
	     	this->print_active_set(active, function, err);
                MSQ_CHKERR(err);
	    });

	    /* check for equilibrium point */
	    if (active->num_active >= 2) {
		MSQ_DEBUG_PRINT(3,"Testing for an equilibrium point \n");
                this->check_equilibrium(&equilibrium_pt, &opt_status, err); 
                       MSQ_CHKERR(err);

		MSQ_DEBUG_ACTION(2,{
		    if (equilibrium_pt) 
			fprintf(stdout,"Optimization Exiting: An equilibrium point \n");
                });
	    }

	    /* record the values */
            current_active_value = active->true_active_value;
	    prev_active_values[iter_count] = active->true_active_value;

	} else {
	    /* decrease the iteration count by one */
	    /* smooth_param->iter_count -= 1; */
	    iter_count -= 1;
	    MSQ_DEBUG_ACTION(2,{
		fprintf(stdout,"Optimization Exiting: No viable directions; equilibrium point \n");
		/* Print the old active set */
		this->print_active_set(active,function,err); MSQ_CHKERR(err);
	    });
	}
      }

      MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Checking the validity of the mesh\n");
	    if (!this->validity_check(err)) fprintf(stdout,"The final mesh is not valid\n");
       MSQ_CHKERR(err);
      });

      MSQ_DEBUG_ACTION(2,{fprintf(stdout,"Number of optimization iterations %d\n",
                            iter_count);});
 
      switch(opt_status) {
	case MSQ_EQUILIBRIUM:
	  MSQ_DEBUG_PRINT(2,"Optimization Termination Opt_Status: Equilibrium\n"); break;
	case MSQ_STEP_TOO_SMALL:
	  MSQ_DEBUG_PRINT(2,"Optimization Termination Opt_Status: Step Too Small\n"); break;
	case MSQ_IMP_TOO_SMALL:
	  MSQ_DEBUG_PRINT(2,"Optimization Termination Opt_Status: Improvement Too Small\n"); break;
	case MSQ_FLAT_NO_IMP:
	  MSQ_DEBUG_PRINT(2,"Optimization Termination Opt_Status: Flat No Improvement\n"); break;
	case MSQ_ZERO_SEARCH:
	  MSQ_DEBUG_PRINT(2,"Optimization Termination Opt_Status: Zero Search\n"); break;
	case MSQ_MAX_ITER_EXCEEDED:
	  MSQ_DEBUG_PRINT(2,"Optimization Termination Opt_Status: Max Iter Exceeded\n"); break;
      }
      
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::step_acceptance"
void NonSmoothSteepestDescent::step_acceptance(PatchData &pd, MsqError &err)
{
//  int        ierr;
//  int        step_accepted;
  int        i;
  int        num_values, num_steps;
  int        valid, step_status;
  int        accept_alpha;
  double     estimated_improvement;
  double     current_improvement = 1E300;
  double     previous_improvement = 1E300;
  double     current_percent_diff = 1E300;
  double     original_point[3];

  num_values = num_function_values;

//  step_accepted = 0;
  step_status = MSQ_STEP_NOT_DONE;
  opt_status = 0;
  num_steps = 0;

  if (alpha < min_step_size) {
      opt_status = MSQ_IMP_TOO_SMALL;
      step_status = MSQ_STEP_DONE;
      MSQ_DEBUG_PRINT(3,"Alpha starts too small, no improvement\n");
  }

  /* save the original function and active set */
  MSQ_COPY_VECTOR(original_point,coords[0],dimension);
  MSQ_COPY_VECTOR(original_function, function, num_values);
  this->copy_active(active, original_active, err); MSQ_CHKERR(err);

  while (step_status == MSQ_STEP_NOT_DONE) {

    num_steps++;  if (num_steps >= 100) step_status = MSQ_STEP_DONE;

    accept_alpha = MSQ_FALSE;

    while (!accept_alpha && alpha>min_step_size) {

      /* make the step */
      for (i=0;i<dimension;i++) {
         coords[0][i] -= alpha*search[i];
      }
        //pd.set_coords_array_element(coords[0],0,err);

      MSQ_DEBUG_ACTION(2,{
         fprintf(stdout,"search direction %f %f \n",search[0],search[1]); 
         fprintf(stdout,"new vertex position %f %f \n",coords[0][0],coords[0][1]); 
      });

      /* assume alpha is acceptable */
      accept_alpha=MSQ_TRUE;

      /* never take a step that makes a valid mesh invalid or worsens the quality */
      valid = validity_check(err); MSQ_CHKERR(err);
      if (valid) valid=improvement_check(err); MSQ_CHKERR(err);
      if (!valid) {
          accept_alpha=MSQ_FALSE;
          for (i=0;i<dimension;i++) {
             coords[0][i] += alpha*search[i];
          }
            //pd.set_coords_array_element(coords[0],0,err);
          alpha = alpha/2;
           MSQ_DEBUG_ACTION(2,{
               fprintf(stdout,"Step not accepted, the new alpha %f\n",alpha); 
          });

          if (alpha < min_step_size) {
 	        opt_status = MSQ_STEP_TOO_SMALL;
                step_status = MSQ_STEP_DONE;
                MSQ_DEBUG_PRINT(2,"Step too small\n");
 	        /* get back the original point, function, and active set */
                MSQ_COPY_VECTOR(coords[0],original_point,dimension);
                  //pd.set_coords_array_element(coords[0],0,err);
	        MSQ_COPY_VECTOR(function,original_function,num_values);
	        this->copy_active(original_active, active, err); 
	  }
       }
    } 
         
    if (valid  && (alpha > min_step_size)) {
      /* compute the new function and active set */
      this->compute_function(&pd, function, err); MSQ_CHKERR(err);
      this->find_active_set(function, active, err); MSQ_CHKERR(err);
	
      /* estimate the minimum improvement by taking this step */
      this->get_min_estimate(&estimated_improvement, err); MSQ_CHKERR(err);
      MSQ_DEBUG_ACTION(3,{
           fprintf(stdout,"The estimated improvement for this step: %f\n",
		   estimated_improvement); 
      });
	
      /* calculate the actual increase */
      current_improvement = active->true_active_value - prev_active_values[iter_count-1];

      MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Actual improvement %f\n",current_improvement);});

      /* calculate the percent difference from estimated increase */
      current_percent_diff = fabs(current_improvement-estimated_improvement)/
	fabs(estimated_improvement);

      /* determine whether to accept a step */
      if ((fabs(previous_improvement) > fabs(current_improvement)) && 
	  (previous_improvement < 0)) {
	/* accept the previous step - it was better */
	     MSQ_DEBUG_PRINT(2,"Accepting the previous step\n");
 
	/* subtract alpha in again (previous step) */
	for (i=0;i<dimension;i++) {
	  coords[0][i] -= alpha*search[i];
	}
            //pd.set_coords_array_element(coords[0],0,err);

	/* does this make an invalid mesh valid? */
        valid = 1;
        valid = validity_check(err); MSQ_CHKERR(err);
        if (valid) valid=improvement_check(err); MSQ_CHKERR(err);

	/* copy test function and active set */
	MSQ_COPY_VECTOR(function,test_function,num_function_values);
	this->copy_active(test_active, active, err); MSQ_CHKERR(err);
 
	opt_status = MSQ_STEP_ACCEPTED;  step_status = MSQ_STEP_DONE;
            
	/* check to see that we're still making good improvements */
	if (fabs(previous_improvement) < min_acceptable_improvement) {
	  opt_status = MSQ_IMP_TOO_SMALL; step_status = MSQ_STEP_DONE;
	  MSQ_DEBUG_PRINT(2,"Optimization Exiting: Improvement too small\n");
	}

      } else if (((fabs(current_improvement) > fabs(estimated_improvement)) ||
		  (current_percent_diff < .1)) && (current_improvement<0)) {
	/* accept this step, exceeded estimated increase or was close */
	opt_status = MSQ_STEP_ACCEPTED;  step_status = MSQ_STEP_DONE;

	/* check to see that we're still making good improvements */
	if (fabs(current_improvement) < min_acceptable_improvement) {
	  MSQ_DEBUG_PRINT(2,"Optimization Exiting: Improvement too small\n");
	  opt_status = MSQ_IMP_TOO_SMALL; step_status = MSQ_STEP_DONE;
	}

      } else if ((current_improvement > 0) && (previous_improvement > 0) &&
		 (fabs(current_improvement) < min_acceptable_improvement) &&
		 (fabs(previous_improvement) < min_acceptable_improvement)) {

	/* we are making no progress, quit */
	opt_status = MSQ_FLAT_NO_IMP; step_status = MSQ_STEP_DONE;
	MSQ_DEBUG_PRINT(2,"Opimization Exiting: Flat no improvement\n");
           
	/* get back the original point, function, and active set */
	MSQ_COPY_VECTOR(coords[0],original_point,dimension);
            //pd.set_coords_array_element(coords[0],0,err);
	MSQ_COPY_VECTOR(function,original_function,num_function_values);
	this->copy_active(original_active, active, err); MSQ_CHKERR(err);

      }
      else
      {
	/* halve alpha and try again */
	/* add out the old step */
	for (i=0;i<dimension;i++) coords[0][i] += alpha*search[i];
            //pd.set_coords_array_element(coords[0],0,err);

	/* halve step size */
	alpha = alpha/2; 
	MSQ_DEBUG_ACTION(3,{fprintf(stdout,"Step not accepted, the new alpha %f\n",alpha); });

	if (alpha < min_step_size)
          {
	  /* get back the original point, function, and active set */
	  MSQ_DEBUG_PRINT(2,"Optimization Exiting: Step too small\n");
	  MSQ_COPY_VECTOR(coords[0],original_point,dimension);
              //pd.set_coords_array_element(coords[0],0,err);
	  MSQ_COPY_VECTOR(function,original_function,num_function_values);
	  this->copy_active(original_active, active, err); MSQ_CHKERR(err);
	  opt_status = MSQ_STEP_TOO_SMALL;  step_status = MSQ_STEP_DONE;
	}
          else
          {
	  MSQ_COPY_VECTOR(test_function, function, num_function_values);
	  this->copy_active(active, test_active, err); MSQ_CHKERR(err);
	  previous_improvement = current_improvement;
	}
      }
    }
  }
  if (current_improvement>0 && opt_status==MSQ_STEP_ACCEPTED) {
    MSQ_DEBUG_ACTION(2,{printf("Accepted a negative step %f \n",current_improvement);});
  }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::get_min_estimate"
void  NonSmoothSteepestDescent::get_min_estimate(double *final_est, MsqError &err)
{
    int    i;
    double est_imp;

    *final_est = -1E300;
    for (i=0;i<active->num_active;i++) {
	est_imp = -alpha*gs[active->active_ind[i]];
        if (est_imp>*final_est) *final_est = est_imp;
    }
    if (*final_est == 0) {
	*final_est = -1E300;
	for (i=0;i<num_function_values;i++) {
	    est_imp = -alpha*gs[i];
	    if ((est_imp>*final_est) && (fabs(est_imp) > MSQ_MACHINE_EPS)) {
		*final_est = est_imp;
	    }
	}
    }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::get_gradient_projections"
void NonSmoothSteepestDescent::get_gradient_projections(MsqError &err)
{
    for (int i=0;i<num_function_values;i++) 
	MSQ_DOT(gs[i],gradient[i],search,dimension);

    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"steepest in get_gradient_projections %d\n",steepest);});
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::compute_alpha"
void NonSmoothSteepestDescent::compute_alpha(MsqError &err)
{
//    int       ierr;
//    int       j;
    int       i;
//    int       ind;
    int       num_values;
    double    steepest_function;
    double    steepest_grad;
    double    alpha_i;
    double    min_positive_value=1E300;

    MSQ_DEBUG_PRINT(2,"In compute alpha\n");

    num_values = num_function_values;
    alpha = 1E300;

    steepest_function = function[steepest];
    steepest_grad = gs[steepest];
    for (i=0;i<num_values;i++)
    {
        /* if it's not active */
      if (i!=steepest)
      {
	  alpha_i = steepest_function - function[i];
	   
	  if (fabs(gs[steepest] - gs[i])>1E-13) {
	     /* compute line intersection */
	     alpha_i = alpha_i/(steepest_grad - gs[i]);
	  } else {
	     /* the lines don't intersect - it's not under consideration*/
	     alpha_i = 0;
	  }
	  if ((alpha_i > min_step_size ) && (fabs(alpha_i) < fabs(alpha))) {
	    alpha = fabs(alpha_i); 
            MSQ_DEBUG_ACTION(3,{ fprintf(stdout,"Setting alpha %d %g\n",i,alpha_i); });
	  }
          if ((alpha_i > 0) && (alpha_i < min_positive_value)) {
            min_positive_value = alpha_i;
          }
       }
    }

    if ((alpha == 1E300) && (min_positive_value != 1E300)) {
      alpha = min_positive_value;
    }

    /* if it never gets set, set it to the default */
    if (alpha == 1E300) {
      alpha = max_alpha;
      MSQ_DEBUG_ACTION(3,{ fprintf(stdout,"Setting alpha to the maximum step length\n"); });
    }

    MSQ_DEBUG_ACTION(3,{ fprintf(stdout,"  The initial step size: %f\n",alpha); });
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::copy_active"
void NonSmoothSteepestDescent::copy_active(ActiveSet *active1, ActiveSet *active2, 
                                          MsqError &err)
{
    if (active1==NULL || active2==NULL)
       err.set_msg("Null memory in copy_active\n");

    active2->num_active = active1->num_active;
    active2->num_equal  = active1->num_equal;
    active2->true_active_value = active1->true_active_value;
    for (int i=0;i<active1->num_active;i++) {
	active2->active_ind[i] = active1->active_ind[i];
    }
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::print_active_set" 
void NonSmoothSteepestDescent::print_active_set(ActiveSet *active_set, 
                                        double *func, MsqError &err)
{
    if (active_set==0) err.set_msg("Null ActiveSet \n");
    if (active_set->num_active == 0) std::cout<< "No active values\n";
    /* print the active set */
    for (int i=0;i<active_set->num_active;i++) {
           fprintf(stdout,"Active value %d:   %f \n",
	               i+1,func[active_set->active_ind[i]]); 
    }
}

#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::init_opt"
void NonSmoothSteepestDescent::init_opt(MsqError &err)
{
    int        i, j;

    MSQ_DEBUG_PRINT(2,"\nInitializing Optimization \n");
    if (num_function_values > 150) {
      err.set_msg("num_values exceeds 150");
    }
    /* for the purposes of initialization will be set to zero after */
    equilibrium_pt = 0;
    step_too_small = 0;
//    step_accepted = 0;
    opt_status = 0;
    iter_count = 0;
    opt_iter_count = 0;
    max_iterations = 100;
    steepest = 0;
    alpha = 0;
    max_alpha = 0;

    MSQ_DEBUG_PRINT(3,"  Initialized Constants \n");
    for (i=0;i<3;i++) {
      search[i] = 0;
      PDG_ind[i] = -1;
      for (j=0;j<3;j++) PDG[i][j] = 0;
    }

    MSQ_DEBUG_PRINT(3,"  Initialized search and PDG \n");
    for (i=0;i<num_function_values;i++) {
       function[i] = 0;
       test_function[i] = 0;
       original_function[i] = 0;
       gs[i] = 0;
       for (j=0;j<3;j++) {
           gradient[i][j] = 0;
       }
    }
    MSQ_DEBUG_PRINT(3,"  Initialized function/gradient \n");
    if (num_function_values > 150) {
      for (i=0;i<150;i++) {
       for (j=0;j<150;j++) G[i][j] = -1;
      }
    } else {
      for (i=0;i<num_function_values;i++) {
       for (j=0;j<num_function_values;j++) G[i][j] = -1;
      }
    }
    MSQ_DEBUG_PRINT(3,"  Initialized G\n");
 
    for (i=0;i<100;i++) prev_active_values[i] = 0;
    MSQ_DEBUG_PRINT(3,"  Initialized prev_active_values\n");
}


#undef __FUNC__
#define __FUNC__ "NonSmoothSteepestDescent::init_max_step_length"
void NonSmoothSteepestDescent::init_max_step_length(MsqError &err)
{
  int i, j, k;
  double max_diff = 0;
  double diff=0;

  MSQ_DEBUG_PRINT(2,"In init_max_step_length\n");

  /* check that the input data is correct */
  if (numElements==0) err.set_msg("Num incident vtx = 0\n");
  if ((dimension!=2) && (dimension!=3)) {
     err.set_msg("Problem dimension is incorrect\n");
  }

  /* find the maximum distance between two incident vertex locations */
  for (i=1;i<numVertices;i++) {
    for (j=i;j<numVertices+1;j++) {
      diff=0;
      for (k=0;k<dimension;k++) {
        diff += (coords[i][k]-coords[j][k])*(coords[i][k]-coords[j][k]);
      }
      if (max_diff < diff) max_diff=diff;
    } 
  }
  max_diff = sqrt(max_diff);
  if (max_diff==0) {
     err.set_msg("Maximum distance between incident vertices = 0\n");
  }
  max_alpha = max_diff/100;

  MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Maximum step is %g\n",max_alpha);});
}

