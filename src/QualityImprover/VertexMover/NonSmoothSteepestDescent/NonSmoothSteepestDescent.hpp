/*!
  \file   NonSmoothSteepestDescent.hpp
  \brief  

  The NonSmoothSteepestDescent Class implements the steepest descent algorythm in
  order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QaulityMetric object.

  \author Thomas Leurent
  \date   2002-06-13
*/

#ifndef Mesquite_NonSmoothSteepestDescent_hpp 
#define Mesquite_NonSmoothSteepestDescent_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "ObjectiveFunction.hpp"

namespace Mesquite
{
#define MSQ_XDIR 0
#define MSQ_YDIR 1
#define MSQ_ZDIR 2
#define MSQ_BIG_POS_NMBR    1E300
#define MSQ_BIG_NEG_NMBR   -1E300
#define MSQ_MAX_OPT_ITER    20

#define MSQ_CCW               1
#define MSQ_CW                0
#define MSQ_NO_EQUIL 101
#define MSQ_CHECK_TOP_DOWN 102
#define MSQ_CHECK_BOTTOM_UP 103
#define MSQ_TWO_PT_PLANE 104
#define MSQ_THREE_PT_PLANE 105
#define MSQ_CHECK_Y_COORD_DIRECTION 106
#define MSQ_CHECK_X_COORD_DIRECTION 107
#define MSQ_CHECK_Z_COORD_DIRECTION 108
#define MSQ_EQUIL 109
#define MSQ_HULL_TEST_ERROR 110

#define MSQ_STEP_ACCEPTED     100
#define MSQ_IMP_TOO_SMALL     101
#define MSQ_FLAT_NO_IMP       102
#define MSQ_STEP_TOO_SMALL    103
#define MSQ_EQUILIBRIUM       104
#define MSQ_ZERO_SEARCH       105
#define MSQ_MAX_ITER_EXCEEDED 106

#define MSQ_STEP_DONE        101
#define MSQ_STEP_NOT_DONE    102

#define MAX_NUM_ELEMENTS 150
#define MAX_FUNC_PER_ELEMENT 6
#define MSQ_MACHINE_EPS     1E-16
#define MSQ_TRUE  1
#define MSQ_FALSE 0
#define MSQ_MAX(a,b) (a > b ? a : b)
#define MSQ_MIN(a,b) (a < b ? a : b)
#define MSQ_LESS_THAN_MACHINE_EPS(x)   ( ((fabs(x)+1.0) > 1.0) ? 0 : 1 )

#define MSQ_DOT(c,a,b,n) {\
  int i99; \
  if (n==2) c = a[0]*b[0] + a[1]*b[1]; \
  else if (n==3) c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];\
  else { \
    for (i99=0;i99<n;i99++) c += a[i99]*b[i99]; \
  } \
}

#define MSQ_NORMALIZE(v,n) {\
    int i99; \
    double mag99; \
    if (n==2){ \
       mag99 = sqrt(v[0]*v[0] + v[1]*v[1]) ; \
       if (mag99 != 0) { \
          v[0] = v[0]/mag99; \
          v[1] = v[1]/mag99; \
       } \
    } else if (n==3) {\
     mag99 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) ; \
     if (mag99 != 0) { \
         v[0] = v[0]/mag99; \
         v[1] = v[1]/mag99; \
         v[2] = v[2]/mag99; \
     } \
   } else { \
     for (i99=0;i99<n;i99++) mag99+=v[i99]+v[i99]; \
     if (mag99 != 0) { \
       for (i99=0;i99<n;i99++) v[i99] = v[i99]/mag99;\
     } \
   }\
}

#define MSQ_COPY_VECTOR(a,b,n) { \
  int i99; \
  if (n==2) { \
     a[0] = b[0];  a[1] = b[1];  \
  } else if (n==3) {\
     a[0] = b[0];  a[1] = b[1];  a[2] = b[2]; \
  } else { \
     for (i99=0;i99<n;i99++) a[i99] = b[i99]; \
  } \
}


  struct ActiveSet
  {
    int num_active;
    int num_equal;
    double true_active_value;
    int active_ind[150];  // need a better way of setting max number of active values
  };

  class NonSmoothSteepestDescent : public VertexMover 
  {
  public:
    NonSmoothSteepestDescent(ObjectiveFunction* of);
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();
    
  private:
    ObjectiveFunction* objFunc;
    
      /* local copy of patch data */
      //    PatchData patch_data;
    int dimension;
    int numVertices;
    int numElements;
    MsqVertex* coords;
    MsqMeshEntity* connectivity;
    
      /* smoothing parameters */
    int max_iterations;
    double conv_eps;
    double active_epsilon;
    double min_acceptable_improvement;
    double min_step_size;
    
      /* optimization data */
    double original_value;
    int iter_count;
    int opt_iter_count;
    int num_function_values;
    double *function;
    double *test_function;
    double *original_function;
    double **gradient;
    int patch_validity;
    int opt_status;
    int equilibrium_pt;
    int step_too_small;
    int step_accepted;
    int max_iter;
    int steepest;
    double search[3];
    double alpha;
    double max_alpha;
    double *gs;
    double *prev_active_values;
    double **G;
    double **PDG;
    int PDG_ind[3];
    int num_LI;
    ActiveSet *active;
    ActiveSet *test_active;
    ActiveSet *original_active;
    double current_active_value;
    
      /* functions */
    void init_opt(MsqError &err);
    void init_max_step_length(MsqError &err);
    
      /* optimize */
    void minmax_opt(PatchData &pd, MsqError &err);
    void step_acceptance(PatchData &pd, MsqError &err); 
    void get_min_estimate(double *final_est, MsqError &err);
    void get_gradient_projections(MsqError &err);
    void compute_alpha(MsqError &err);
    void copy_active(ActiveSet *active1, ActiveSet *active2, MsqError &err);
    
      /* function/gradient/active set computations */
    void compute_function(PatchData *pd, double *function, MsqError &err);
    double** compute_gradient(PatchData *pd, MsqError &err);
    void find_active_set(double *function, ActiveSet *active_set, MsqError &err);
    void print_active_set(ActiveSet *active_set, double *func, MsqError &err);
    
      /* checking validity/improvement */
    int improvement_check(MsqError &err);
    int validity_check(MsqError &err);
    
      /* checking equilibrium routines */
    void check_equilibrium(int *equil, int *opt_status, MsqError &err);
    int convex_hull_test(double **vec, int num_vec, MsqError &err);
    int check_vector_dots(double **vec, int num_vec, double *normal, MsqError &err);
    void find_plane_normal(double pt1[3], double pt2[3], double pt3[3], 
                           double *cross, MsqError &err);
    void find_plane_points(int dir1, int dir2, double **vec, int num_vec, double *pt1,
                           double *pt2, double*pt3, int *opt_status, MsqError &err);
    
      /* from the matrix file */
    void form_grammian(double **vec, MsqError &err);
    void form_PD_grammian(MsqError &err);
    void singular_test(int n, double **A, int *singular, MsqError &err);
    void condition3x3(double **A, double *cond, MsqError &err);
    void solve2x2(double a11, double a12, double a21, double a22, 
                  double b1, double b2, double **x,MsqError &err);
    void form_reduced_matrix(double ***P, MsqError &err);
    
      /* search direction */
    void search_direction(PatchData &pd, MsqError &err);
    void search_edges_faces(double **dir, MsqError &err);
    void get_active_directions(double **gradient,
                               double ***dir, MsqError &err);
  };
  
}

#endif
