/*!
  \file   MeanRatioQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \author Todd Munson
  \author Thomas Leurent

  \date   2002-06-9
*/
#include <vector>
#include "MeanRatioQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"
#include "MsqTimer.hpp"

using namespace Mesquite;
using std::cout;
using std::endl;

#undef __FUNC__
#define __FUNC__ "MeanRatioQualityMetric::MeanRatioQualityMetric"

MeanRatioQualityMetric::MeanRatioQualityMetric()
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
  set_negate_flag(1);
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  set_name("Mean Ratio");
}

/*****************************************************************************/
/* Not all compilers substitute out constants (especially the square root).  */
/* Therefore, they are substituted out manually.  The values below were      */
/* calculated on a solaris machine using long doubles. I believe they are    */
/* accurate.                                                                 */
/*****************************************************************************/

#define isqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0)*/
#define tisqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0)*/
#define isqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0)*/
#define tisqrt6  1.22474487139158915752946794252e+00        /*  3.0/sqrt(6.0)*/

#define a2       5.00000000000000000000000000000e-01        /*  1.0/2.0      */
#define b2      -1.00000000000000000000000000000e-00        /* -1.0/1.0      */
#define b2m1    -2.00000000000000000000000000000e-00        /* -2.0/1.0      */

#define a3       3.33333333333333333333333333333e-01        /*  1.0/3.0      */
#define b3      -6.66666666666666666666666666667e-01        /* -2.0/3.0      */
#define b3m1    -1.66666666666666666666666666667e-00        /* -5.0/3.0      */

/*****************************************************************************/
/* This set of functions reference triangular elements to an equilateral     */
/* triangle.  The input are the coordinates in the following order:          */
/*      [x1 x2 x3 y1 y2 y3]                                                  */
/* A unit normal is also required in order to calculate the value.           */
/* A zero return value indicates success, while a nonzero value indicates    */
/* failure.                                                                  */
/*****************************************************************************/

/*****************************************************************************/
/* Function evaluation -- requires 43 flops.                                 */
/*****************************************************************************/
inline bool m_fcn_2e(double &obj, const Vector3D x[3], const Vector3D &n)
{
  double matr[9], f;
  double g;

  /* Calculate M = [A*inv(W) n] */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = (2.0*x[2][0] - x[1][0] - x[0][0])*isqrt3;
  matr[2] = n[0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = (2.0*x[2][1] - x[1][1] - x[0][1])*isqrt3;
  matr[5] = n[1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = (2.0*x[2][2] - x[1][2] - x[0][2])*isqrt3;
  matr[8] = n[2];

  /* Calculate det(M). */
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[3]*(matr[2]*matr[7] - matr[1]*matr[8]) +
      matr[6]*(matr[1]*matr[5] - matr[2]*matr[4]);
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] +
      matr[3]*matr[3] + matr[4]*matr[4] +
      matr[6]*matr[6] + matr[7]*matr[7];

  /* Calculate objective function. */
  obj = a2 * f * pow(g, b2);
  return true;
}

/*****************************************************************************/
/* Gradient requires 85 flops.                                               */
/*****************************************************************************/
inline bool g_fcn_2e(double &obj, Vector3D g_obj[3], 
                     const Vector3D x[3], const Vector3D &n)
{
  double matr[9], f;
  double adj_m[9], g;		// adj_m[2,5,8] not used
  double loc1, loc2, loc3, loc4;

  /* Calculate M = [A*inv(W) n] */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = (2.0*x[2][0] - x[1][0] - x[0][0])*isqrt3;
  matr[2] = n[0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = (2.0*x[2][1] - x[1][1] - x[0][1])*isqrt3;
  matr[5] = n[1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = (2.0*x[2][2] - x[1][2] - x[0][2])*isqrt3;
  matr[8] = n[2];

  /* Calculate det([n M]). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[2]*matr[7] - matr[1]*matr[8];
  loc3 = matr[1]*matr[5] - matr[2]*matr[4];
  g = matr[0]*loc1 + matr[3]*loc2 + matr[6]*loc3;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] +
      matr[3]*matr[3] + matr[4]*matr[4] +
      matr[6]*matr[6] + matr[7]*matr[7];

  /* Calculate objective function. */
  loc4 = a2 * pow(g, b2);
  obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc4;
  g = b2 * obj / g;

  adj_m[0] = matr[0]*f + loc1*g;
  adj_m[3] = matr[3]*f + loc2*g;
  adj_m[6] = matr[6]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[3]*g;
  loc3 = matr[6]*g;

  adj_m[1] = matr[1]*f + loc3*matr[5] - loc2*matr[8];
  adj_m[4] = matr[4]*f + loc1*matr[8] - loc3*matr[2];
  adj_m[7] = matr[7]*f + loc2*matr[2] - loc1*matr[5];

  loc1 = isqrt3*adj_m[1];
  g_obj[0][0] = -adj_m[0] - loc1;
  g_obj[1][0] =  adj_m[0] - loc1;
  g_obj[2][0] = 2.0*loc1;

  loc1 = isqrt3*adj_m[4];
  g_obj[0][1] = -adj_m[3] - loc1;
  g_obj[1][1] =  adj_m[3] - loc1;
  g_obj[2][1] = 2.0*loc1;

  loc1 = isqrt3*adj_m[7];
  g_obj[0][2] = -adj_m[6] - loc1;
  g_obj[1][2] =  adj_m[6] - loc1;
  g_obj[2][2] = 2.0*loc1;
  return true;
}

/*****************************************************************************/
/* Hessian requires 289 flops.                                               */
/*****************************************************************************/
inline bool h_fcn_2e(double &obj, Vector3D g_obj[3], Matrix3D h_obj[6],
                     const Vector3D x[3], const Vector3D &n)
{
  double matr[9], f;
  double adj_m[9], g;		// adj_m[2,5,8] not used
  double dg[9];			// dg[2,5,8] not used
  double loc0, loc1, loc2, loc3, loc4;
  double A[12], J_A[6], J_B[9], J_C[9];  // only 2x2 corners used

  /* Calculate M = [A*inv(W) n] */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = (2.0*x[2][0] - x[1][0] - x[0][0])*isqrt3;
  matr[2] = n[0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = (2.0*x[2][1] - x[1][1] - x[0][1])*isqrt3;
  matr[5] = n[1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = (2.0*x[2][2] - x[1][2] - x[0][2])*isqrt3;
  matr[8] = n[2];

  /* Calculate det([n M]). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  g = matr[0]*dg[0] + matr[3]*dg[3] + matr[6]*dg[6];
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] +
      matr[3]*matr[3] + matr[4]*matr[4] +
      matr[6]*matr[6] + matr[7]*matr[7];

  loc4 = g;

  /* Calculate objective function. */
  loc1 = a2 * pow(g, b2);
  obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc1;
  g = b2 * obj / g;

  dg[3] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];

  adj_m[0] = matr[0]*f + dg[0]*g;
  adj_m[1] = matr[1]*f + dg[1]*g;
  adj_m[3] = matr[3]*f + dg[3]*g;
  adj_m[4] = matr[4]*f + dg[4]*g;
  adj_m[6] = matr[6]*f + dg[6]*g;
  adj_m[7] = matr[7]*f + dg[7]*g;

  loc1 = isqrt3*adj_m[1];
  g_obj[0][0] = -adj_m[0] - loc1;
  g_obj[1][0] =  adj_m[0] - loc1;
  g_obj[2][0] = 2.0*loc1;

  loc1 = isqrt3*adj_m[4];
  g_obj[0][1] = -adj_m[3] - loc1;
  g_obj[1][1] =  adj_m[3] - loc1;
  g_obj[2][1] = 2.0*loc1;

  loc1 = isqrt3*adj_m[7];
  g_obj[0][2] = -adj_m[6] - loc1;
  g_obj[1][2] =  adj_m[6] - loc1;
  g_obj[2][2] = 2.0*loc1;

  loc0 = g;
  loc1 = f;
  f = f*b2/loc4;
  g = g*b2m1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = dg[0]*f;
  loc4 = dg[0]*g + loc2;

  J_A[0] = loc1 + dg[0]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] + loc4*dg[1];
  J_B[0] = loc3*matr[3] + loc4*dg[3];
  J_B[1] = loc3*matr[4] + loc4*dg[4];
  J_C[0] = loc3*matr[6] + loc4*dg[6];
  J_C[1] = loc3*matr[7] + loc4*dg[7];

  loc2 = matr[1]*f;
  loc3 = dg[1]*f;
  loc4 = dg[1]*g + loc2;

  J_A[3] = loc1 + dg[1]*(loc2 + loc4);
  J_B[3] = loc3*matr[3] + loc4*dg[3];
  J_B[4] = loc3*matr[4] + loc4*dg[4];
  J_C[3] = loc3*matr[6] + loc4*dg[6];
  J_C[4] = loc3*matr[7] + loc4*dg[7];

  /* First diagonal block */
  loc2 = isqrt3*J_A[1];
  A[0] = -J_A[0] - loc2;
  A[1] =  J_A[0] - loc2;

  loc2 = isqrt3*J_A[3];
  A[4] = -J_A[1] - loc2;
  A[5] =  J_A[1] - loc2;
  A[6] = 2.0*loc2;

  loc2 = isqrt3*A[4];
  h_obj[0][0][0] = -A[0] - loc2;
  h_obj[1][0][0] =  A[0] - loc2;
  h_obj[2][0][0] = 2.0*loc2;

  loc2 = isqrt3*A[5];
  h_obj[3][0][0] = A[1] - loc2;
  h_obj[4][0][0] = 2.0*loc2;

  h_obj[5][0][0] = tisqrt3*A[6];

  /* First off-diagonal block */
  loc2 = matr[8]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = isqrt3*J_B[3];
  A[0] = -J_B[0] - loc2;
  A[1] =  J_B[0] - loc2;
  A[2] = 2.0*loc2;

  loc2 = isqrt3*J_B[4];
  A[4] = -J_B[1] - loc2;
  A[5] =  J_B[1] - loc2;
  A[6] = 2.0*loc2;

  loc2 = isqrt3*A[4];
  h_obj[0][0][1] = -A[0] - loc2;
  h_obj[1][0][1] =  A[0] - loc2;
  h_obj[2][0][1] = 2.0*loc2;

  loc2 = isqrt3*A[5];
  h_obj[1][1][0] = -A[1] - loc2;
  h_obj[3][0][1] =  A[1] - loc2;
  h_obj[4][0][1] = 2.0*loc2;

  loc2 = isqrt3*A[6];
  h_obj[2][1][0] = -A[2] - loc2;
  h_obj[4][1][0] =  A[2] - loc2;
  h_obj[5][0][1] = 2.0*loc2;

  /* Second off-diagonal block */
  loc2 = matr[5]*loc0;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  loc2 = isqrt3*J_C[3];
  A[0] = -J_C[0] - loc2;
  A[1] =  J_C[0] - loc2;
  A[2] = 2.0*loc2;

  loc2 = isqrt3*J_C[4];
  A[4] = -J_C[1] - loc2;
  A[5] =  J_C[1] - loc2;
  A[6] = 2.0*loc2;

  loc2 = isqrt3*A[4];
  h_obj[0][0][2] = -A[0] - loc2;
  h_obj[1][0][2] =  A[0] - loc2;
  h_obj[2][0][2] = 2.0*loc2;

  loc2 = isqrt3*A[5];
  h_obj[1][2][0] = -A[1] - loc2;
  h_obj[3][0][2] =  A[1] - loc2;
  h_obj[4][0][2] = 2.0*loc2;

  loc2 = isqrt3*A[6];
  h_obj[2][2][0] = -A[2] - loc2;
  h_obj[4][2][0] =  A[2] - loc2;
  h_obj[5][0][2] = 2.0*loc2;

  /* Second block of rows */
  loc2 = matr[3]*f;
  loc3 = dg[3]*f;
  loc4 = dg[3]*g + loc2;

  J_A[0] = loc1 + dg[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[4] + loc4*dg[4];
  J_B[0] = loc3*matr[6] + loc4*dg[6];
  J_B[1] = loc3*matr[7] + loc4*dg[7];

  loc2 = matr[4]*f;
  loc3 = dg[4]*f;
  loc4 = dg[4]*g + loc2;

  J_A[3] = loc1 + dg[4]*(loc2 + loc4);
  J_B[3] = loc3*matr[6] + loc4*dg[6];
  J_B[4] = loc3*matr[7] + loc4*dg[7];

  /* Second diagonal block */
  loc2 = isqrt3*J_A[1];
  A[0] = -J_A[0] - loc2;
  A[1] =  J_A[0] - loc2;

  loc2 = isqrt3*J_A[3];
  A[4] = -J_A[1] - loc2;
  A[5] =  J_A[1] - loc2;
  A[6] = 2.0*loc2;

  loc2 = isqrt3*A[4];
  h_obj[0][1][1] = -A[0] - loc2;
  h_obj[1][1][1] =  A[0] - loc2;
  h_obj[2][1][1] = 2.0*loc2;

  loc2 = isqrt3*A[5];
  h_obj[3][1][1] = A[1] - loc2;
  h_obj[4][1][1] = 2.0*loc2;

  h_obj[5][1][1] = tisqrt3*A[6];

  /* Third off-diagonal block */
  loc2 = matr[2]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = isqrt3*J_B[3];
  A[0] = -J_B[0] - loc2;
  A[1] =  J_B[0] - loc2;
  A[2] = 2.0*loc2;

  loc2 = isqrt3*J_B[4];
  A[4] = -J_B[1] - loc2;
  A[5] =  J_B[1] - loc2;
  A[6] = 2.0*loc2;

  loc2 = isqrt3*A[4];
  h_obj[0][1][2] = -A[0] - loc2;
  h_obj[1][1][2] =  A[0] - loc2;
  h_obj[2][1][2] = 2.0*loc2;

  loc2 = isqrt3*A[5];
  h_obj[1][2][1] = -A[1] - loc2;
  h_obj[3][1][2] =  A[1] - loc2;
  h_obj[4][1][2] = 2.0*loc2;

  loc2 = isqrt3*A[6];
  h_obj[2][2][1] = -A[2] - loc2;
  h_obj[4][2][1] =  A[2] - loc2;
  h_obj[5][1][2] = 2.0*loc2;

  /* Third block of rows */
  loc2 = matr[6]*f;
  loc3 = dg[6]*f;
  loc4 = dg[6]*g + loc2;

  J_A[0] = loc1 + dg[6]*(loc2 + loc4);
  J_A[1] = loc3*matr[7] + loc4*dg[7];

  loc2 = matr[7]*f;
  loc4 = dg[7]*g + loc2;

  J_A[3] = loc1 + dg[7]*(loc2 + loc4);

  /* Third diagonal block */
  loc2 = isqrt3*J_A[1];
  A[0] = -J_A[0] - loc2;
  A[1] =  J_A[0] - loc2;

  loc2 = isqrt3*J_A[3];
  A[4] = -J_A[1] - loc2;
  A[5] =  J_A[1] - loc2;
  A[6] = 2.0*loc2;

  loc2 = isqrt3*A[4];
  h_obj[0][2][2] = -A[0] - loc2;
  h_obj[1][2][2] =  A[0] - loc2;
  h_obj[2][2][2] = 2.0*loc2;

  loc2 = isqrt3*A[5];
  h_obj[3][2][2] = A[1] - loc2;
  h_obj[4][2][2] = 2.0*loc2;

  h_obj[5][2][2] = tisqrt3*A[6];

  // completes diagonal blocks.
  h_obj[0].fill_lower_triangle();
  h_obj[3].fill_lower_triangle();
  h_obj[5].fill_lower_triangle();
  return true;
}

/*****************************************************************************/
/* This set of functions reference triangular elements to an isoscles        */
/* triangle.  The input are the coordinates in the following order:          */
/*      [x1 x2 x3 y1 y2 y3]                                                  */
/* A unit normal is also required in order to calculate the value.           */
/* A zero return value indicates success, while a nonzero value indicates    */
/* failure.                                                                  */
/*****************************************************************************/

/*****************************************************************************/
/* Function evaluation -- requires 34 flops.                                 */
/*****************************************************************************/
inline bool m_fcn_2i(double &obj, const Vector3D x[3], const Vector3D &n)
{
  double matr[9], f;
  double g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = x[2][0] - x[0][0];
  matr[2] = n[0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = x[2][1] - x[0][1];
  matr[5] = n[1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = x[2][2] - x[0][2];
  matr[8] = n[2];

  /* Calculate det(M). */
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[3]*(matr[2]*matr[7] - matr[1]*matr[8]) +
      matr[6]*(matr[1]*matr[5] - matr[2]*matr[4]);
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] +
      matr[3]*matr[3] + matr[4]*matr[4] +
      matr[6]*matr[6] + matr[7]*matr[7];

  /* Calculate objective function. */
  obj = a2 * f * pow(g, b2);
  return true;
}

/*****************************************************************************/
/* Gradient requires 67 flops.                                               */
/*****************************************************************************/
inline bool g_fcn_2i(double &obj, Vector3D g_obj[3], 
                     const Vector3D x[3], const Vector3D &n)
{
  double matr[9], f;
  double adj_m[9], g;            // adj_m[2,5,8] not used
  double loc1, loc2, loc3, loc4;

  /* Calculate M = [A*inv(W) n] */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = x[2][0] - x[0][0];
  matr[2] = n[0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = x[2][1] - x[0][1];
  matr[5] = n[1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = x[2][2] - x[0][2];
  matr[8] = n[2];

  /* Calculate det([n M]). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[2]*matr[7] - matr[1]*matr[8];
  loc3 = matr[1]*matr[5] - matr[2]*matr[4];
  g = matr[0]*loc1 + matr[3]*loc2 + matr[6]*loc3;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] +
      matr[3]*matr[3] + matr[4]*matr[4] +
      matr[6]*matr[6] + matr[7]*matr[7];

  /* Calculate objective function. */
  loc4 = a2 * pow(g, b2);
  obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc4;
  g = b2 * obj / g;

  adj_m[0] = matr[0]*f + loc1*g;
  adj_m[3] = matr[3]*f + loc2*g;
  adj_m[6] = matr[6]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[3]*g;
  loc3 = matr[6]*g;

  adj_m[1] = matr[1]*f + loc3*matr[5] - loc2*matr[8];
  adj_m[4] = matr[4]*f + loc1*matr[8] - loc3*matr[2];
  adj_m[7] = matr[7]*f + loc2*matr[2] - loc1*matr[5];

  g_obj[0][0] = -adj_m[0] - adj_m[1];
  g_obj[1][0] =  adj_m[0];
  g_obj[2][0] =  adj_m[1];

  g_obj[0][1] = -adj_m[3] - adj_m[4];
  g_obj[1][1] =  adj_m[3];
  g_obj[2][1] =  adj_m[4];

  g_obj[0][2] = -adj_m[6] - adj_m[7];
  g_obj[1][2] =  adj_m[6];
  g_obj[2][2] =  adj_m[7];
  return true;
}

/*****************************************************************************/
/* Hessian requires 190 flops.                                               */
/*****************************************************************************/
inline bool h_fcn_2i(double &obj, Vector3D g_obj[3], Matrix3D h_obj[6],
                     const Vector3D x[3], const Vector3D &n)
{
  double matr[9], f;
  double adj_m[9], g;		// adj_m[2,5,8] not used
  double dg[9];			// dg[2,5,8] not used
  double loc0, loc1, loc2, loc3, loc4;
  double A[12], J_A[6], J_B[9], J_C[9];  // only 2x2 corners used

  /* Calculate M = [A*inv(W) n] */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = x[2][0] - x[0][0];
  matr[2] = n[0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = x[2][1] - x[0][1];
  matr[5] = n[1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = x[2][2] - x[0][2];
  matr[8] = n[2];

  /* Calculate det([n M]). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  g = matr[0]*dg[0] + matr[3]*dg[3] + matr[6]*dg[6];
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] +
      matr[3]*matr[3] + matr[4]*matr[4] +
      matr[6]*matr[6] + matr[7]*matr[7];

  loc4 = g;

  /* Calculate objective function. */
  loc1 = a2 * pow(g, b2);
  obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc1;
  g = b2 * obj / g;

  dg[3] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];

  adj_m[0] = matr[0]*f + dg[0]*g;
  adj_m[1] = matr[1]*f + dg[1]*g;
  adj_m[3] = matr[3]*f + dg[3]*g;
  adj_m[4] = matr[4]*f + dg[4]*g;
  adj_m[6] = matr[6]*f + dg[6]*g;
  adj_m[7] = matr[7]*f + dg[7]*g;

  g_obj[0][0] = -adj_m[0] - adj_m[1];
  g_obj[1][0] =  adj_m[0];
  g_obj[2][0] =  adj_m[1];

  g_obj[0][1] = -adj_m[3] - adj_m[4];
  g_obj[1][1] =  adj_m[3];
  g_obj[2][1] =  adj_m[4];

  g_obj[0][2] = -adj_m[6] - adj_m[7];
  g_obj[1][2] =  adj_m[6];
  g_obj[2][2] =  adj_m[7];

  loc0 = g;
  loc1 = f;
  f = f*b2/loc4;
  g = g*b2m1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = dg[0]*f;
  loc4 = dg[0]*g + loc2;

  J_A[0] = loc1 + dg[0]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] + loc4*dg[1];
  J_B[0] = loc3*matr[3] + loc4*dg[3];
  J_B[1] = loc3*matr[4] + loc4*dg[4];
  J_C[0] = loc3*matr[6] + loc4*dg[6];
  J_C[1] = loc3*matr[7] + loc4*dg[7];

  loc2 = matr[1]*f;
  loc3 = dg[1]*f;
  loc4 = dg[1]*g + loc2;

  J_A[3] = loc1 + dg[1]*(loc2 + loc4);
  J_B[3] = loc3*matr[3] + loc4*dg[3];
  J_B[4] = loc3*matr[4] + loc4*dg[4];
  J_C[3] = loc3*matr[6] + loc4*dg[6];
  J_C[4] = loc3*matr[7] + loc4*dg[7];

  /* First diagonal block */
  A[0] = -J_A[0] - J_A[1];
  A[1] =  J_A[0];

  A[4] = -J_A[1] - J_A[3];
  A[5] =  J_A[1];
  A[6] =  J_A[3];

  h_obj[0][0][0] = -A[0] - A[4];
  h_obj[1][0][0] =  A[0];
  h_obj[2][0][0] =  A[4];

  h_obj[3][0][0] = A[1];
  h_obj[4][0][0] = A[5];

  h_obj[5][0][0] = A[6];

  /* First off-diagonal block */
  loc2 = matr[8]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  A[0] = -J_B[0] - J_B[3];
  A[1] =  J_B[0];
  A[2] =  J_B[3];

  A[4] = -J_B[1] - J_B[4];
  A[5] =  J_B[1];
  A[6] =  J_B[4];

  h_obj[0][0][1] = -A[0] - A[4];
  h_obj[1][0][1] =  A[0];
  h_obj[2][0][1] =  A[4];

  h_obj[1][1][0] = -A[1] - A[5];
  h_obj[3][0][1] =  A[1];
  h_obj[4][0][1] =  A[5];

  h_obj[2][1][0] = -A[2] - A[6];
  h_obj[4][1][0] =  A[2];
  h_obj[5][0][1] =  A[6];

  /* Second off-diagonal block */
  loc2 = matr[5]*loc0;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  A[0] = -J_C[0] - J_C[3];
  A[1] =  J_C[0];
  A[2] =  J_C[3];

  A[4] = -J_C[1] - J_C[4];
  A[5] =  J_C[1];
  A[6] =  J_C[4];

  h_obj[0][0][2] = -A[0] - A[4];
  h_obj[1][0][2] =  A[0];
  h_obj[2][0][2] =  A[4];

  h_obj[1][2][0] = -A[1] - A[5];
  h_obj[3][0][2] =  A[1];
  h_obj[4][0][2] =  A[5];

  h_obj[2][2][0] = -A[2] - A[6];
  h_obj[4][2][0] =  A[2];
  h_obj[5][0][2] =  A[6];

  /* Second block of rows */
  loc2 = matr[3]*f;
  loc3 = dg[3]*f;
  loc4 = dg[3]*g + loc2;

  J_A[0] = loc1 + dg[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[4] + loc4*dg[4];
  J_B[0] = loc3*matr[6] + loc4*dg[6];
  J_B[1] = loc3*matr[7] + loc4*dg[7];

  loc2 = matr[4]*f;
  loc3 = dg[4]*f;
  loc4 = dg[4]*g + loc2;

  J_A[3] = loc1 + dg[4]*(loc2 + loc4);
  J_B[3] = loc3*matr[6] + loc4*dg[6];
  J_B[4] = loc3*matr[7] + loc4*dg[7];

  /* Second diagonal block */
  A[0] = -J_A[0] - J_A[1];
  A[1] =  J_A[0];

  A[4] = -J_A[1] - J_A[3];
  A[5] =  J_A[1];
  A[6] =  J_A[3];

  h_obj[0][1][1] = -A[0] - A[4];
  h_obj[1][1][1] =  A[0];
  h_obj[2][1][1] =  A[4];

  h_obj[3][1][1] = A[1];
  h_obj[4][1][1] = A[5];

  h_obj[5][1][1] = A[6];

  /* Third off-diagonal block */
  loc2 = matr[2]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  A[0] = -J_B[0] - J_B[3];
  A[1] =  J_B[0];
  A[2] =  J_B[3];

  A[4] = -J_B[1] - J_B[4];
  A[5] =  J_B[1];
  A[6] =  J_B[4];

  h_obj[0][1][2] = -A[0] - A[4];
  h_obj[1][1][2] =  A[0];
  h_obj[2][1][2] =  A[4];

  h_obj[1][2][1] = -A[1] - A[5];
  h_obj[3][1][2] =  A[1];
  h_obj[4][1][2] =  A[5];

  h_obj[2][2][1] = -A[2] - A[6];
  h_obj[4][2][1] =  A[2];
  h_obj[5][1][2] =  A[6];

  /* Third block of rows */
  loc2 = matr[6]*f;
  loc3 = dg[6]*f;
  loc4 = dg[6]*g + loc2;

  J_A[0] = loc1 + dg[6]*(loc2 + loc4);
  J_A[1] = loc3*matr[7] + loc4*dg[7];

  loc2 = matr[7]*f;
  loc4 = dg[7]*g + loc2;

  J_A[3] = loc1 + dg[7]*(loc2 + loc4);

  /* Third diagonal block */
  A[0] = -J_A[0] - J_A[1];
  A[1] =  J_A[0];

  A[4] = -J_A[1] - J_A[3];
  A[5] =  J_A[1];
  A[6] =  J_A[3];

  h_obj[0][2][2] = -A[0] - A[4];
  h_obj[1][2][2] =  A[0];
  h_obj[2][2][2] =  A[4];

  h_obj[3][2][2] = A[1];
  h_obj[4][2][2] = A[5];

  h_obj[5][2][2] = A[6];

  // completes diagonal blocks.
  h_obj[0].fill_lower_triangle();
  h_obj[3].fill_lower_triangle();
  h_obj[5].fill_lower_triangle();
  return true;
}

/*****************************************************************************/
/* The following set of functions reference tetrahedral elements to a        */
/* regular tetrahedron.  They are used when assessing the quality of a       */
/* tetrahedral element.  A zero return value indicates success, while        */
/* a nonzero value indicates failure.                                        */
/*****************************************************************************/

/*****************************************************************************/
/* Function evaluation -- requires 61 flops.                                 */
/*****************************************************************************/
inline bool m_fcn_3e(double &obj, const Vector3D x[4])
{
  double matr[9], f;
  double g;

  /* Calculate M = A*inv(W). */
  f       = x[1][0] + x[0][0];
  matr[0] = x[1][0] - x[0][0];
  matr[1] = (2.0*x[2][0] - f)*isqrt3;
  matr[2] = (3.0*x[3][0] - x[2][0] - f)*isqrt6;

  f       = x[1][1] + x[0][1];
  matr[3] = x[1][1] - x[0][1];
  matr[4] = (2.0*x[2][1] - f)*isqrt3;
  matr[5] = (3.0*x[3][1] - x[2][1] - f)*isqrt6;

  f       = x[1][2] + x[0][2];
  matr[6] = x[1][2] - x[0][2];
  matr[7] = (2.0*x[2][2] - f)*isqrt3;
  matr[8] = (3.0*x[3][2] - x[2][2] - f)*isqrt6;

  /* Calculate det(M). */
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  obj = a3 * f * pow(g, b3);
  return true;
}

/*****************************************************************************/
/* Optimal gradient calculation courtesy of Paul Hovland (at least we        */
/* think it is optimal).  The original code provided was modified to         */
/* reduce the number of flops and intermediate variables, and improve the    */
/* locality of reference.                                                    */
/*                                                                           */
/* This function requires 130 flops.                                         */
/*****************************************************************************/

inline bool g_fcn_3e(double &obj, Vector3D g_obj[4], const Vector3D x[4])
{
  double matr[9], f;
  double adj_m[9], g;
  double loc1, loc2, loc3, loc4;

  /* Calculate M = A*inv(W). */
  f       = x[1][0] + x[0][0];
  matr[0] = x[1][0] - x[0][0];
  matr[1] = (2.0*x[2][0] - f)*isqrt3;
  matr[2] = (3.0*x[3][0] - x[2][0] - f)*isqrt6;

  f       = x[1][1] + x[0][1];
  matr[3] = x[1][1] - x[0][1];
  matr[4] = (2.0*x[2][1] - f)*isqrt3;
  matr[5] = (3.0*x[3][1] - x[2][1] - f)*isqrt6;

  f       = x[1][2] + x[0][2];
  matr[6] = x[1][2] - x[0][2];
  matr[7] = (2.0*x[2][2] - f)*isqrt3;
  matr[8] = (3.0*x[3][2] - x[2][2] - f)*isqrt6;

  /* Calculate det(M). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[5]*matr[6] - matr[3]*matr[8];
  loc3 = matr[3]*matr[7] - matr[4]*matr[6];
  g = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  loc4 = a3 * pow(g, b3);
  obj  = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc4;
  g = b3 * obj / g;

  adj_m[0] = matr[0]*f + loc1*g;
  adj_m[1] = matr[1]*f + loc2*g;
  adj_m[2] = matr[2]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[1]*g;
  loc3 = matr[2]*g;

  adj_m[3] = matr[3]*f + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = matr[4]*f + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = matr[5]*f + loc2*matr[6] - loc1*matr[7];

  adj_m[6] = matr[6]*f + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = matr[7]*f + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = matr[8]*f + loc1*matr[4] - loc2*matr[3];

  loc1 = isqrt3*adj_m[1];
  loc2 = isqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0][0] = -adj_m[0] - loc3;
  g_obj[1][0] = adj_m[0] - loc3;
  g_obj[2][0] = 2.0*loc1 - loc2;
  g_obj[3][0] = 3.0*loc2;

  loc1 = isqrt3*adj_m[4];
  loc2 = isqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[0][1] = -adj_m[3] - loc3;
  g_obj[1][1] = adj_m[3] - loc3;
  g_obj[2][1] = 2.0*loc1 - loc2;
  g_obj[3][1] = 3.0*loc2;

  loc1 = isqrt3*adj_m[7];
  loc2 = isqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[0][2] = -adj_m[6] - loc3;
  g_obj[1][2] = adj_m[6] - loc3;
  g_obj[2][2] = 2.0*loc1 - loc2;
  g_obj[3][2] = 3.0*loc2;
  return true;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1-d3) each contain 10 elements, while the  */
/* off-diagonal elements (b1-b3) each contain 16 elements.                   */
/*                                                                           */
/* This function requires 598 flops.                                         */
/*****************************************************************************/
/* The form of the function, gradient, and Hessian is the following:         */
/*   o(x) = a * f(A(x)) * pow(g(A(x)), b)                                    */
/* where A(x) is the matrix generated from:                                  */
/*           [x1-x0 x2-x0 x3-x0]                                             */
/*    A(x) = [y1-y0 y2-y0 y3-y0] * inv(W)                                    */
/*           [z1-z0 z2-z0 z3-z0]                                             */
/* and f() is the squared Frobenius norm of A(x), and g() is the determinant */
/* of A(x).                                                                  */
/*                                                                           */
/* The gradient is calculated as follows:                                    */
/*   alpha := a*pow(g(A(x)),b)                                               */
/*   beta  := a*b*f(A(x))*pow(g(A(x)),b-1)                                   */
/*                                                                           */
/*                                                                           */
/*   do/dx = (alpha * (df/dA) + beta * (dg/dA)) (dA/dx)                      */
/*                                                                           */
/*   Note: this is the optimal ordering for the gradient vector.             */
/*   Distributing (dA/dx) would result in two matrix vector products as      */
/*   opposed to the single matrix vector product in the above formulation.   */
/*                                                                           */
/*   (df/dA)_i = 2*A_i                                                       */
/*   (dg/dA)_i = A_j*A_k - A_l*A_m for some {j,k,l,m}                        */
/*                                                                           */
/*   d^2o/dx^2 = (dA/dx)' * ((d alpha/dA) * (df/dA) +                        */
/*                           (d  beta/dA) * (dg/dA)                          */
/*                                  alpha * (d^2f/dA^2)                      */
/*                                   beta * (d^2g/dA^2)) * (dA/dx)           */
/*                                                                           */
/*   Note: since A(x) is a linear function, there are no terms involving     */
/*   d^2A/dx^2 since this matrix is zero.                                    */
/*                                                                           */
/*   gamma := a*b*pow(g(A(x)),b-1)                                           */
/*   delta := a*b*(b-1)*f(A(x))*pow(g(A(x)),b-2)                             */
/*                                                                           */
/*   d^2o/dx^2 = (dA/dx)' * (gamma*((dg/dA)'*(df/dA) + (df/dA)'*(dg/dA)) +   */
/*                           delta* (dg/dA)'*(dg/dA) +                       */
/*                           alpha*(d^2f/dA^2) +                             */
/*                            beta*(d^2g/dA^2)) * (dA/dx)                    */
/*                                                                           */
/*   Note: (df/dA) and (dg/dA) are row vectors and we only calculate the     */
/*   upper triangular part of the inner matrix.                              */
/*                                                                           */
/*   For regular tetrahedral elements, we have the following:                */
/*                                                                           */
/*           [-1         1        0         0         ]                      */
/*       M = [-sqrt(3)  -sqrt(3)  2*sqrt(3) 0         ]                      */
/*           [-sqrt(6)  -sqrt(6)  -sqrt(6)  3*sqrt(6) ]                      */
/*                                                                           */
/*           [M 0 0]                                                         */
/*   dA/dx = [0 M 0]                                                         */
/*           [0 0 M]                                                         */
/*                                                                           */
/*   I belive the above is close to optimal for the calculation of the       */
/*   Hessian.  Distributing the (dA/dx) results in larger vector which are   */
/*   detrimental when forming the outer product.  The way the method is      */
/*   written, we only calculate a 9x9 symmetric matrix in the outer product. */
/*                                                                           */
/*   In two dimensions, the inner matrix computed has a nice structure and   */
/*   we can eliminate some of the computation in the inner product.  This    */
/*   does not appear to be the case in more than two dimensions.             */
/*****************************************************************************/

inline bool h_fcn_3e(double &obj, Vector3D g_obj[4], Matrix3D h_obj[10], 
		     const Vector3D x[4])
{
  double matr[9], f;
  double adj_m[9], g;
  double dg[9], loc0, loc1, loc2, loc3, loc4;
  double A[12], J_A[6], J_B[9], J_C[9];

  /* Calculate M = A*inv(W). */
  f       = x[1][0] + x[0][0];
  matr[0] = x[1][0] - x[0][0];
  matr[1] = (2.0*x[2][0] - f)*isqrt3;
  matr[2] = (3.0*x[3][0] - x[2][0] - f)*isqrt6;

  f       = x[1][1] + x[0][1];
  matr[3] = x[1][1] - x[0][1];
  matr[4] = (2.0*x[2][1] - f)*isqrt3;
  matr[5] = (3.0*x[3][1] - x[2][1] - f)*isqrt6;

  f       = x[1][2] + x[0][2];
  matr[6] = x[1][2] - x[0][2];
  matr[7] = (2.0*x[2][2] - f)*isqrt3;
  matr[8] = (3.0*x[3][2] - x[2][2] - f)*isqrt6;

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  g = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  loc4 = g;

  /* Calculate objective function. */
  loc1 = a3 * pow(g, b3);
  obj  = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc1;
  g = b3 * obj / g; 

  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  adj_m[0] = matr[0]*f + dg[0]*g;
  adj_m[1] = matr[1]*f + dg[1]*g;
  adj_m[2] = matr[2]*f + dg[2]*g;
  adj_m[3] = matr[3]*f + dg[3]*g;
  adj_m[4] = matr[4]*f + dg[4]*g;
  adj_m[5] = matr[5]*f + dg[5]*g;
  adj_m[6] = matr[6]*f + dg[6]*g;
  adj_m[7] = matr[7]*f + dg[7]*g;
  adj_m[8] = matr[8]*f + dg[8]*g;

  loc1 = isqrt3*adj_m[1];
  loc2 = isqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0][0] = -adj_m[0] - loc3;
  g_obj[1][0] = adj_m[0] - loc3;
  g_obj[2][0] = 2.0*loc1 - loc2;
  g_obj[3][0] = 3.0*loc2;

  loc1 = isqrt3*adj_m[4];
  loc2 = isqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[0][1] = -adj_m[3] - loc3;
  g_obj[1][1] = adj_m[3] - loc3;
  g_obj[2][1] = 2.0*loc1 - loc2;
  g_obj[3][1] = 3.0*loc2;

  loc1 = isqrt3*adj_m[7];
  loc2 = isqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[0][2] = -adj_m[6] - loc3;
  g_obj[1][2] = adj_m[6] - loc3;
  g_obj[2][2] = 2.0*loc1 - loc2;
  g_obj[3][2] = 3.0*loc2;

  loc0 = g;
  loc1 = f;
  f = f*b3/loc4;
  g = g*b3m1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = dg[0]*f;
  loc4 = dg[0]*g + loc2;

  J_A[0] = loc1 + dg[0]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] + loc4*dg[1];
  J_A[2] = loc3*matr[2] + loc4*dg[2];
  J_B[0] = loc3*matr[3] + loc4*dg[3];
  J_B[1] = loc3*matr[4] + loc4*dg[4];
  J_B[2] = loc3*matr[5] + loc4*dg[5];
  J_C[0] = loc3*matr[6] + loc4*dg[6];
  J_C[1] = loc3*matr[7] + loc4*dg[7];
  J_C[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[1]*f;
  loc3 = dg[1]*f;
  loc4 = dg[1]*g + loc2;

  J_A[3] = loc1 + dg[1]*(loc2 + loc4);
  J_A[4] = loc3*matr[2] + loc4*dg[2];
  J_B[3] = loc3*matr[3] + loc4*dg[3];
  J_B[4] = loc3*matr[4] + loc4*dg[4];
  J_B[5] = loc3*matr[5] + loc4*dg[5];
  J_C[3] = loc3*matr[6] + loc4*dg[6];
  J_C[4] = loc3*matr[7] + loc4*dg[7];
  J_C[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[2]*f;
  loc3 = dg[2]*f;
  loc4 = dg[2]*g + loc2;

  J_A[5] = loc1 + dg[2]*(loc2 + loc4);
  J_B[6] = loc3*matr[3] + loc4*dg[3];
  J_B[7] = loc3*matr[4] + loc4*dg[4];
  J_B[8] = loc3*matr[5] + loc4*dg[5];
  J_C[6] = loc3*matr[6] + loc4*dg[6];
  J_C[7] = loc3*matr[7] + loc4*dg[7];
  J_C[8] = loc3*matr[8] + loc4*dg[8];

  /* First diagonal block */
  loc2 = isqrt3*J_A[1];
  loc3 = isqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = isqrt3*J_A[3];
  loc3 = isqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = isqrt3*J_A[4];
  loc3 = isqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = isqrt3*A[4];
  loc3 = isqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0][0][0] = -A[0] - loc4;
  h_obj[1][0][0] =  A[0] - loc4;
  h_obj[2][0][0] = 2.0*loc2 - loc3;
  h_obj[3][0][0] = 3.0*loc3;

  loc2 = isqrt3*A[5];
  loc3 = isqrt6*A[9];

  h_obj[4][0][0] = A[1] - loc2 - loc3;
  h_obj[5][0][0] = 2.0*loc2 - loc3;
  h_obj[6][0][0] = 3.0*loc3;

  loc3 = isqrt6*A[10];
  h_obj[7][0][0] = tisqrt3*A[6] - loc3;
  h_obj[8][0][0] = 3.0*loc3;

  h_obj[9][0][0] = tisqrt6*A[11];

  /* First off-diagonal block */
  loc2 = matr[8]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[7]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[6]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = isqrt3*J_B[3];
  loc3 = isqrt6*J_B[6];
  loc4 = loc2 + loc3;

  A[0] = -J_B[0] - loc4;
  A[1] =  J_B[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = isqrt3*J_B[4];
  loc3 = isqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;
  A[5] =  J_B[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = isqrt3*J_B[5];
  loc3 = isqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;
  A[9] =  J_B[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = isqrt3*A[4];
  loc3 = isqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0][0][1] = -A[0] - loc4;
  h_obj[1][0][1] =  A[0] - loc4;
  h_obj[2][0][1] = 2.0*loc2 - loc3;
  h_obj[3][0][1] = 3.0*loc3;

  loc2 = isqrt3*A[5];
  loc3 = isqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[1][1][0] = -A[1] - loc4;
  h_obj[4][0][1] =  A[1] - loc4;
  h_obj[5][0][1] = 2.0*loc2 - loc3;
  h_obj[6][0][1] = 3.0*loc3;

  loc2 = isqrt3*A[6];
  loc3 = isqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[2][1][0] = -A[2] - loc4;
  h_obj[5][1][0] =  A[2] - loc4;
  h_obj[7][0][1] = 2.0*loc2 - loc3;
  h_obj[8][0][1] = 3.0*loc3;

  loc2 = isqrt3*A[7];
  loc3 = isqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[3][1][0] = -A[3] - loc4;
  h_obj[6][1][0] =  A[3] - loc4;
  h_obj[8][1][0] = 2.0*loc2 - loc3;
  h_obj[9][0][1] = 3.0*loc3;

  /* Second off-diagonal block */
  loc2 = matr[5]*loc0;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  loc2 = matr[4]*loc0;
  J_C[2] += loc2;
  J_C[6] -= loc2;

  loc2 = matr[3]*loc0;
  J_C[5] -= loc2;
  J_C[7] += loc2;

  loc2 = isqrt3*J_C[3];
  loc3 = isqrt6*J_C[6];
  loc4 = loc2 + loc3;

  A[0] = -J_C[0] - loc4;
  A[1] =  J_C[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = isqrt3*J_C[4];
  loc3 = isqrt6*J_C[7];
  loc4 = loc2 + loc3;

  A[4] = -J_C[1] - loc4;
  A[5] =  J_C[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = isqrt3*J_C[5];
  loc3 = isqrt6*J_C[8];
  loc4 = loc2 + loc3;

  A[8] = -J_C[2] - loc4;
  A[9] =  J_C[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = isqrt3*A[4];
  loc3 = isqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0][0][2] = -A[0] - loc4;
  h_obj[1][0][2] =  A[0] - loc4;
  h_obj[2][0][2] = 2.0*loc2 - loc3;
  h_obj[3][0][2] = 3.0*loc3;

  loc2 = isqrt3*A[5];
  loc3 = isqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[1][2][0] = -A[1] - loc4;
  h_obj[4][0][2] =  A[1] - loc4;
  h_obj[5][0][2] = 2.0*loc2 - loc3;
  h_obj[6][0][2] = 3.0*loc3;

  loc2 = isqrt3*A[6];
  loc3 = isqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[2][2][0] = -A[2] - loc4;
  h_obj[5][2][0] =  A[2] - loc4;
  h_obj[7][0][2] = 2.0*loc2 - loc3;
  h_obj[8][0][2] = 3.0*loc3;

  loc2 = isqrt3*A[7];
  loc3 = isqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[3][2][0] = -A[3] - loc4;
  h_obj[6][2][0] =  A[3] - loc4;
  h_obj[8][2][0] = 2.0*loc2 - loc3;
  h_obj[9][0][2] = 3.0*loc3;

  /* Second block of rows */
  loc2 = matr[3]*f;
  loc3 = dg[3]*f;
  loc4 = dg[3]*g + loc2;

  J_A[0] = loc1 + dg[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[4] + loc4*dg[4];
  J_A[2] = loc3*matr[5] + loc4*dg[5];
  J_B[0] = loc3*matr[6] + loc4*dg[6];
  J_B[1] = loc3*matr[7] + loc4*dg[7];
  J_B[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[4]*f;
  loc3 = dg[4]*f;
  loc4 = dg[4]*g + loc2;

  J_A[3] = loc1 + dg[4]*(loc2 + loc4);
  J_A[4] = loc3*matr[5] + loc4*dg[5];
  J_B[3] = loc3*matr[6] + loc4*dg[6];
  J_B[4] = loc3*matr[7] + loc4*dg[7];
  J_B[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[5]*f;
  loc3 = dg[5]*f;
  loc4 = dg[5]*g + loc2;

  J_A[5] = loc1 + dg[5]*(loc2 + loc4);
  J_B[6] = loc3*matr[6] + loc4*dg[6];
  J_B[7] = loc3*matr[7] + loc4*dg[7];
  J_B[8] = loc3*matr[8] + loc4*dg[8];

  /* Second diagonal block */
  loc2 = isqrt3*J_A[1];
  loc3 = isqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = isqrt3*J_A[3];
  loc3 = isqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = isqrt3*J_A[4];
  loc3 = isqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = isqrt3*A[4];
  loc3 = isqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0][1][1] = -A[0] - loc4;
  h_obj[1][1][1] =  A[0] - loc4;
  h_obj[2][1][1] = 2.0*loc2 - loc3;
  h_obj[3][1][1] = 3.0*loc3;

  loc2 = isqrt3*A[5];
  loc3 = isqrt6*A[9];

  h_obj[4][1][1] = A[1] - loc2 - loc3;
  h_obj[5][1][1] = 2.0*loc2 - loc3;
  h_obj[6][1][1] = 3.0*loc3;

  loc3 = isqrt6*A[10];
  h_obj[7][1][1] = tisqrt3*A[6] - loc3;
  h_obj[8][1][1] = 3.0*loc3;

  h_obj[9][1][1] = tisqrt6*A[11];

  /* Third off-diagonal block */
  loc2 = matr[2]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[1]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[0]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = isqrt3*J_B[3];
  loc3 = isqrt6*J_B[6];
  loc4 = loc2 + loc3;

  A[0] = -J_B[0] - loc4;
  A[1] =  J_B[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = isqrt3*J_B[4];
  loc3 = isqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;
  A[5] =  J_B[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = isqrt3*J_B[5];
  loc3 = isqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;
  A[9] =  J_B[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = isqrt3*A[4];
  loc3 = isqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0][1][2] = -A[0] - loc4;
  h_obj[1][1][2] =  A[0] - loc4;
  h_obj[2][1][2] = 2.0*loc2 - loc3;
  h_obj[3][1][2] = 3.0*loc3;

  loc2 = isqrt3*A[5];
  loc3 = isqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[1][2][1] = -A[1] - loc4;
  h_obj[4][1][2] =  A[1] - loc4;
  h_obj[5][1][2] = 2.0*loc2 - loc3;
  h_obj[6][1][2] = 3.0*loc3;

  loc2 = isqrt3*A[6];
  loc3 = isqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[2][2][1] = -A[2] - loc4;
  h_obj[5][2][1] =  A[2] - loc4;
  h_obj[7][1][2] = 2.0*loc2 - loc3;
  h_obj[8][1][2] = 3.0*loc3;

  loc2 = isqrt3*A[7];
  loc3 = isqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[3][2][1] = -A[3] - loc4;
  h_obj[6][2][1] =  A[3] - loc4;
  h_obj[8][2][1] = 2.0*loc2 - loc3;
  h_obj[9][1][2] = 3.0*loc3;

  /* Third block of rows */
  loc2 = matr[6]*f;
  loc3 = dg[6]*f;
  loc4 = dg[6]*g + loc2;

  J_A[0] = loc1 + dg[6]*(loc2 + loc4);
  J_A[1] = loc3*matr[7] + loc4*dg[7];
  J_A[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[7]*f;
  loc3 = dg[7]*f;
  loc4 = dg[7]*g + loc2;

  J_A[3] = loc1 + dg[7]*(loc2 + loc4);
  J_A[4] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[8]*f;
  loc4 = dg[8]*g + loc2;

  J_A[5] = loc1 + dg[8]*(loc2 + loc4);

  /* Third diagonal block */
  loc2 = isqrt3*J_A[1];
  loc3 = isqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = isqrt3*J_A[3];
  loc3 = isqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = isqrt3*J_A[4];
  loc3 = isqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = isqrt3*A[4];
  loc3 = isqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0][2][2] = -A[0] - loc4;
  h_obj[1][2][2] =  A[0] - loc4;
  h_obj[2][2][2] = 2.0*loc2 - loc3;
  h_obj[3][2][2] = 3.0*loc3;

  loc2 = isqrt3*A[5];
  loc3 = isqrt6*A[9];

  h_obj[4][2][2] = A[1] - loc2 - loc3;
  h_obj[5][2][2] = 2.0*loc2 - loc3;
  h_obj[6][2][2] = 3.0*loc3;

  loc3 = isqrt6*A[10];
  h_obj[7][2][2] = tisqrt3*A[6] - loc3;
  h_obj[8][2][2] = 3.0*loc3;

  h_obj[9][2][2] = tisqrt6*A[11];

  // completes diagonal blocks.
  h_obj[0].fill_lower_triangle();
  h_obj[4].fill_lower_triangle();
  h_obj[7].fill_lower_triangle();
  h_obj[9].fill_lower_triangle();
  
  return true;
}

/*****************************************************************************/
/* The following set of functions reference tetrahedral elements to a        */
/* right tetrahedron.  They are used when assessing the quality of a         */
/* hexahedral element.  A zero return value indicates success, while         */
/* a nonzero value indicates failure.                                        */
/*****************************************************************************/

/*****************************************************************************/
/* Function evaluation -- requires 43 flops.                                 */
/*****************************************************************************/

inline bool m_fcn_3i(double &obj, const Vector3D x[4])
{
  double matr[9], f;
  double g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = x[2][0] - x[0][0];
  matr[2] = x[3][0] - x[0][0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = x[2][1] - x[0][1];
  matr[5] = x[3][1] - x[0][1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = x[2][2] - x[0][2];
  matr[8] = x[3][2] - x[0][2];

  /* Calculate det(M). */
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  obj = a3 * f * pow(g, b3);
  return true;
}

/*****************************************************************************/
/* Optimal gradient calculation courtesy of Paul Hovland (at least we        */
/* think it is optimal).  The original code provided was modified to         */
/* reduce the number of flops and intermediate variables, and improve the    */
/* locality of reference.                                                    */
/*                                                                           */
/* This function requires 94 flops.                                          */
/*****************************************************************************/

inline bool g_fcn_3i(double &obj, Vector3D g_obj[4], const Vector3D x[12])
{
  double matr[9], f;
  double adj_m[9], g;
  double loc1, loc2, loc3, loc4;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = x[2][0] - x[0][0];
  matr[2] = x[3][0] - x[0][0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = x[2][1] - x[0][1];
  matr[5] = x[3][1] - x[0][1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = x[2][2] - x[0][2];
  matr[8] = x[3][2] - x[0][2];

  /* Calculate det(M). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[5]*matr[6] - matr[3]*matr[8];
  loc3 = matr[3]*matr[7] - matr[4]*matr[6];
  g = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];
 
  /* Calculate objective function. */
  loc4 = a3 * pow(g, b3);
  obj  = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc4;
  g = b3*obj/g; 

  adj_m[0] = matr[0]*f + loc1*g;
  adj_m[1] = matr[1]*f + loc2*g;
  adj_m[2] = matr[2]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[1]*g;
  loc3 = matr[2]*g;

  adj_m[3] = matr[3]*f + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = matr[4]*f + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = matr[5]*f + loc2*matr[6] - loc1*matr[7];

  adj_m[6] = matr[6]*f + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = matr[7]*f + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = matr[8]*f + loc1*matr[4] - loc2*matr[3];

  g_obj[0][0] = -adj_m[0] - adj_m[1] - adj_m[2];
  g_obj[1][0] =  adj_m[0];
  g_obj[2][0] =  adj_m[1];
  g_obj[3][0] =  adj_m[2];

  g_obj[0][1] = -adj_m[3] - adj_m[4] - adj_m[5];
  g_obj[1][1] =  adj_m[3];
  g_obj[2][1] =  adj_m[4];
  g_obj[3][1] =  adj_m[5];

  g_obj[0][2] = -adj_m[6] - adj_m[7] - adj_m[8];
  g_obj[1][2] =  adj_m[6];
  g_obj[2][2] =  adj_m[7];
  g_obj[3][2] =  adj_m[8];
  return true;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1-d3) each contain 10 elements, while the  */
/* off-diagonal elements (b1-b3) each contain 16 elements.                   */
/*                                                                           */
/* The code requires 364 flops.                                              */
/*****************************************************************************/

inline bool h_fcn_3i(double &obj, Vector3D g_obj[4], Matrix3D h_obj[10],
		     const Vector3D x[4])
{
  double matr[9], f;
  double adj_m[9], g;
  double dg[9], loc0, loc1, loc2, loc3, loc4;
  double A[3], J_A[6], J_B[9], J_C[9];

  /* Calculate M = A*inv(W). */
  matr[0] = x[1][0] - x[0][0];
  matr[1] = x[2][0] - x[0][0];
  matr[2] = x[3][0] - x[0][0];

  matr[3] = x[1][1] - x[0][1];
  matr[4] = x[2][1] - x[0][1];
  matr[5] = x[3][1] - x[0][1];

  matr[6] = x[1][2] - x[0][2];
  matr[7] = x[2][2] - x[0][2];
  matr[8] = x[3][2] - x[0][2];

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  g = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  loc4 = g;

  /* Calculate objective function. */
  loc1 = a3 * pow(g, b3);
  obj  = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc1;
  g = b3*obj/g; 

  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  adj_m[0] = matr[0]*f + dg[0]*g;
  adj_m[1] = matr[1]*f + dg[1]*g;
  adj_m[2] = matr[2]*f + dg[2]*g;
  adj_m[3] = matr[3]*f + dg[3]*g;
  adj_m[4] = matr[4]*f + dg[4]*g;
  adj_m[5] = matr[5]*f + dg[5]*g;
  adj_m[6] = matr[6]*f + dg[6]*g;
  adj_m[7] = matr[7]*f + dg[7]*g;
  adj_m[8] = matr[8]*f + dg[8]*g;

  g_obj[0][0] = -adj_m[0] - adj_m[1] - adj_m[2];
  g_obj[1][0] =  adj_m[0];
  g_obj[2][0] =  adj_m[1];
  g_obj[3][0] =  adj_m[2];

  g_obj[0][1] = -adj_m[3] - adj_m[4] - adj_m[5];
  g_obj[1][1] =  adj_m[3];
  g_obj[2][1] =  adj_m[4];
  g_obj[3][1] =  adj_m[5];

  g_obj[0][2] = -adj_m[6] - adj_m[7] - adj_m[8];
  g_obj[1][2] =  adj_m[6];
  g_obj[2][2] =  adj_m[7];
  g_obj[3][2] =  adj_m[8];

  loc0 = g;
  loc1 = f;
  f = f*b3/loc4;
  g = g*b3m1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = dg[0]*f;
  loc4 = dg[0]*g + loc2;

  J_A[0] = loc1 + dg[0]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] + loc4*dg[1];
  J_A[2] = loc3*matr[2] + loc4*dg[2];
  J_B[0] = loc3*matr[3] + loc4*dg[3];
  J_B[1] = loc3*matr[4] + loc4*dg[4];
  J_B[2] = loc3*matr[5] + loc4*dg[5];
  J_C[0] = loc3*matr[6] + loc4*dg[6];
  J_C[1] = loc3*matr[7] + loc4*dg[7];
  J_C[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[1]*f;
  loc3 = dg[1]*f;
  loc4 = dg[1]*g + loc2;

  J_A[3] = loc1 + dg[1]*(loc2 + loc4);
  J_A[4] = loc3*matr[2] + loc4*dg[2];
  J_B[3] = loc3*matr[3] + loc4*dg[3];
  J_B[4] = loc3*matr[4] + loc4*dg[4];
  J_B[5] = loc3*matr[5] + loc4*dg[5];
  J_C[3] = loc3*matr[6] + loc4*dg[6];
  J_C[4] = loc3*matr[7] + loc4*dg[7];
  J_C[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[2]*f;
  loc3 = dg[2]*f;
  loc4 = dg[2]*g + loc2;

  J_A[5] = loc1 + dg[2]*(loc2 + loc4);
  J_B[6] = loc3*matr[3] + loc4*dg[3];
  J_B[7] = loc3*matr[4] + loc4*dg[4];
  J_B[8] = loc3*matr[5] + loc4*dg[5];
  J_C[6] = loc3*matr[6] + loc4*dg[6];
  J_C[7] = loc3*matr[7] + loc4*dg[7];
  J_C[8] = loc3*matr[8] + loc4*dg[8];

  /* First diagonal block */
  A[0] = -J_A[0] - J_A[1] - J_A[2];
  A[1] = -J_A[1] - J_A[3] - J_A[4];
  A[2] = -J_A[2] - J_A[4] - J_A[5];

  h_obj[0][0][0] = -A[0] - A[1] - A[2];
  h_obj[1][0][0] =  A[0];
  h_obj[2][0][0] =  A[1];
  h_obj[3][0][0] =  A[2];

  h_obj[4][0][0] = J_A[0];
  h_obj[5][0][0] = J_A[1];
  h_obj[6][0][0] = J_A[2];

  h_obj[7][0][0] = J_A[3];
  h_obj[8][0][0] = J_A[4];

  h_obj[9][0][0] = J_A[5];

  /* First off-diagonal block */
  loc2 = matr[8]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[7]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[6]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  A[0] = -J_B[0] - J_B[3] - J_B[6];
  A[1] = -J_B[1] - J_B[4] - J_B[7];
  A[2] = -J_B[2] - J_B[5] - J_B[8];

  h_obj[0][0][1] = -A[0] - A[1] - A[2];
  h_obj[1][0][1] =  A[0];
  h_obj[2][0][1] =  A[1];
  h_obj[3][0][1] =  A[2];

  h_obj[1][1][0] = -J_B[0] - J_B[1] - J_B[2];
  h_obj[4][0][1] =  J_B[0];
  h_obj[5][0][1] =  J_B[1];
  h_obj[6][0][1] =  J_B[2];

  h_obj[2][1][0] = -J_B[3] - J_B[4] - J_B[5];
  h_obj[5][1][0] =  J_B[3];
  h_obj[7][0][1] =  J_B[4];
  h_obj[8][0][1] =  J_B[5];

  h_obj[3][1][0] = -J_B[6] - J_B[7] - J_B[8];
  h_obj[6][1][0] =  J_B[6];
  h_obj[8][1][0] =  J_B[7];
  h_obj[9][0][1] =  J_B[8];

  /* Second off-diagonal block */
  loc2 = matr[5]*loc0;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  loc2 = matr[4]*loc0;
  J_C[2] += loc2;
  J_C[6] -= loc2;

  loc2 = matr[3]*loc0;
  J_C[5] -= loc2;
  J_C[7] += loc2;

  A[0] = -J_C[0] - J_C[3] - J_C[6];
  A[1] = -J_C[1] - J_C[4] - J_C[7];
  A[2] = -J_C[2] - J_C[5] - J_C[8];

  h_obj[0][0][2] = -A[0] - A[1] - A[2];
  h_obj[1][0][2] =  A[0];
  h_obj[2][0][2] =  A[1];
  h_obj[3][0][2] =  A[2];

  h_obj[1][2][0] = -J_C[0] - J_C[1] - J_C[2];
  h_obj[4][0][2] =  J_C[0];
  h_obj[5][0][2] =  J_C[1];
  h_obj[6][0][2] =  J_C[2];

  h_obj[2][2][0] = -J_C[3] - J_C[4] - J_C[5];
  h_obj[5][2][0] =  J_C[3];
  h_obj[7][0][2] =  J_C[4];
  h_obj[8][0][2] =  J_C[5];

  h_obj[3][2][0] = -J_C[6] - J_C[7] - J_C[8];
  h_obj[6][2][0] =  J_C[6];
  h_obj[8][2][0] =  J_C[7];
  h_obj[9][0][2] =  J_C[8];

  /* Second block of rows */
  loc2 = matr[3]*f;
  loc3 = dg[3]*f;
  loc4 = dg[3]*g + loc2;

  J_A[0] = loc1 + dg[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[4] + loc4*dg[4];
  J_A[2] = loc3*matr[5] + loc4*dg[5];
  J_B[0] = loc3*matr[6] + loc4*dg[6];
  J_B[1] = loc3*matr[7] + loc4*dg[7];
  J_B[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[4]*f;
  loc3 = dg[4]*f;
  loc4 = dg[4]*g + loc2;

  J_A[3] = loc1 + dg[4]*(loc2 + loc4);
  J_A[4] = loc3*matr[5] + loc4*dg[5];
  J_B[3] = loc3*matr[6] + loc4*dg[6];
  J_B[4] = loc3*matr[7] + loc4*dg[7];
  J_B[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[5]*f;
  loc3 = dg[5]*f;
  loc4 = dg[5]*g + loc2;

  J_A[5] = loc1 + dg[5]*(loc2 + loc4);
  J_B[6] = loc3*matr[6] + loc4*dg[6];
  J_B[7] = loc3*matr[7] + loc4*dg[7];
  J_B[8] = loc3*matr[8] + loc4*dg[8];

  /* Second diagonal block */
  A[0] = -J_A[0] - J_A[1] - J_A[2];
  A[1] = -J_A[1] - J_A[3] - J_A[4];
  A[2] = -J_A[2] - J_A[4] - J_A[5];

  h_obj[0][1][1] = -A[0] - A[1] - A[2];
  h_obj[1][1][1] =  A[0];
  h_obj[2][1][1] =  A[1];
  h_obj[3][1][1] =  A[2];

  h_obj[4][1][1] = J_A[0];
  h_obj[5][1][1] = J_A[1];
  h_obj[6][1][1] = J_A[2];

  h_obj[7][1][1] = J_A[3];
  h_obj[8][1][1] = J_A[4];

  h_obj[9][1][1] = J_A[5];

  /* Third off-diagonal block */
  loc2 = matr[2]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[1]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[0]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  A[0] = -J_B[0] - J_B[3] - J_B[6];
  A[1] = -J_B[1] - J_B[4] - J_B[7];
  A[2] = -J_B[2] - J_B[5] - J_B[8];

  h_obj[0][1][2] = -A[0] - A[1] - A[2];
  h_obj[1][1][2] =  A[0];
  h_obj[2][1][2] =  A[1];
  h_obj[3][1][2] =  A[2];

  h_obj[1][2][1] = -J_B[0] - J_B[1] - J_B[2];
  h_obj[4][1][2] =  J_B[0];
  h_obj[5][1][2] =  J_B[1];
  h_obj[6][1][2] =  J_B[2];

  h_obj[2][2][1] = -J_B[3] - J_B[4] - J_B[5];
  h_obj[5][2][1] =  J_B[3];
  h_obj[7][1][2] =  J_B[4];
  h_obj[8][1][2] =  J_B[5];

  h_obj[3][2][1] = -J_B[6] - J_B[7] - J_B[8];
  h_obj[6][2][1] =  J_B[6];
  h_obj[8][2][1] =  J_B[7];
  h_obj[9][1][2] =  J_B[8];

  /* Third block of rows */
  loc2 = matr[6]*f;
  loc3 = dg[6]*f;
  loc4 = dg[6]*g + loc2;

  J_A[0] = loc1 + dg[6]*(loc2 + loc4);
  J_A[1] = loc3*matr[7] + loc4*dg[7];
  J_A[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[7]*f;
  loc3 = dg[7]*f;
  loc4 = dg[7]*g + loc2;

  J_A[3] = loc1 + dg[7]*(loc2 + loc4);
  J_A[4] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[8]*f;
  loc4 = dg[8]*g + loc2;

  J_A[5] = loc1 + dg[8]*(loc2 + loc4);

  /* Third diagonal block */
  A[0] = -J_A[0] - J_A[1] - J_A[2];
  A[1] = -J_A[1] - J_A[3] - J_A[4];
  A[2] = -J_A[2] - J_A[4] - J_A[5];

  h_obj[0][2][2] = -A[0] - A[1] - A[2];
  h_obj[1][2][2] =  A[0];
  h_obj[2][2][2] =  A[1];
  h_obj[3][2][2] =  A[2];

  h_obj[4][2][2] = J_A[0];
  h_obj[5][2][2] = J_A[1];
  h_obj[6][2][2] = J_A[2];

  h_obj[7][2][2] = J_A[3];
  h_obj[8][2][2] = J_A[4];

  h_obj[9][2][2] = J_A[5];

  // completes diagonal blocks.
  h_obj[0].fill_lower_triangle();
  h_obj[4].fill_lower_triangle();
  h_obj[7].fill_lower_triangle();
  h_obj[9].fill_lower_triangle();

  return true;
}


#undef __FUNC__
#define __FUNC__ "MeanRatioQualityMetric::evaluate_element" 
bool MeanRatioQualityMetric::evaluate_element(PatchData &pd,
                                              MsqMeshEntity *e,
                                              double &m,
                                              MsqError &err)
{
  EntityTopology topo = e->get_element_type();

  MsqVertex *vertices = pd.get_vertex_array(err);
  const size_t *v_i = e->get_vertex_index_array();

  Vector3D coords[4];		// Coordinates for the (decomposed) elements
  Vector3D n;			// Surface normal for 2D objects
  double metrics[8];

  int locs_hex[8][4] = {{0, 1, 3, 4},	// Hex element descriptions
		        {1, 2, 0, 5},
		        {2, 3, 1, 6},
		        {3, 0, 2, 7},
		        {4, 7, 5, 0},
		        {5, 4, 6, 1},
		        {6, 5, 7, 2},
		        {7, 6, 4, 3}};
  int i;

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_CHKERR(err);
    n = n / n.length();		// Need unit normal
    coords[0] = vertices[v_i[0]];
    coords[1] = vertices[v_i[1]];
    coords[2] = vertices[v_i[2]];
    if (!m_fcn_2e(m, coords, n)) return false;
    break;
    
  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_CHKERR(err);
    for (i = 0; i < 4; ++i) {
      n = n / n.length();	// Need unit normal
      coords[0] = vertices[v_i[locs_hex[i][0]]];
      coords[1] = vertices[v_i[locs_hex[i][1]]];
      coords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!m_fcn_2i(metrics[i], coords, n)) return false;
    }
    m = average_metrics(metrics, 4, err);
    break;

  case TETRAHEDRON:
    coords[0] = vertices[v_i[0]];
    coords[1] = vertices[v_i[1]];
    coords[2] = vertices[v_i[2]];
    coords[3] = vertices[v_i[3]];
    if (!m_fcn_3e(m, coords)) return false;
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      coords[0] = vertices[v_i[locs_hex[i][0]]];
      coords[1] = vertices[v_i[locs_hex[i][1]]];
      coords[2] = vertices[v_i[locs_hex[i][2]]];
      coords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!m_fcn_3i(metrics[i], coords)) return false;
    }
    m = average_metrics(metrics, 8, err);
    break;

  default:
    break;
  } // end switch over element type
  return true;
}

#undef __FUNC__
#define __FUNC__ "MeanRatioQualityMetric::compute_element_analytical_gradient" 
bool MeanRatioQualityMetric::compute_element_analytical_gradient(PatchData &pd,
								 MsqMeshEntity *e,
								 MsqVertex *v[], 
								 Vector3D g[],
								 int nv, 
								 double &m,
                                                                 MsqError &err)
{
//  FUNCTION_TIMER_START(__FUNC__);
  
  EntityTopology topo = e->get_element_type();

  if (((topo == QUADRILATERAL) || (topo == HEXAHEDRON)) && 
      ((avgMethod == MINIMUM) || (avgMethod == MAXIMUM))) {
    std::cout << "Minimum and maximum not continuously differentiable." << std::endl;
    std::cout << "Element of subdifferential will be returned." << std::endl;
  }

  MsqVertex *vertices = pd.get_vertex_array(err);
  const size_t *v_i = e->get_vertex_index_array();

  Vector3D coords[4];		// Vertex coordinates for the (decomposed) elements
  Vector3D n;			// Surface normal for 2D objects

  Vector3D gradients[32];	// Gradient of metric with respect to the coords
  Vector3D grad[8];		// Accumulated gradients (composed merit function)
  double   metrics[8];		// Metric values for the (decomposed) elements
  double   nm, t=0;

  int locs_hex[8][4] = {{0, 1, 3, 4},	// Hex element descriptions
                        {1, 2, 0, 5},
                        {2, 3, 1, 6},
                        {3, 0, 2, 7},
                        {4, 7, 5, 0},
                        {5, 4, 6, 1},
                        {6, 5, 7, 2},
                        {7, 6, 4, 3}};
  int i, j;

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_CHKERR(err);
    n = n / n.length();		// Need unit normal
    coords[0] = vertices[v_i[0]];
    coords[1] = vertices[v_i[1]];
    coords[2] = vertices[v_i[2]];
    if (!g_fcn_2e(m, grad, coords, n)) return false;

    // This is not very efficient, but is one way to select correct gradients.
    // For gradients, info is returned only for free vertices, in the 
    // order of v[].
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < nv; ++j) {
        if (vertices + v_i[i] == v[j]) {
          g[j] = grad[i];
        }
      }
    }
    break;

  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_CHKERR(err);
    for (i = 0; i < 4; ++i) {
      grad[i] = 0.0;

      n = n / n.length();	// Need unit normal
      coords[0] = vertices[v_i[locs_hex[i][0]]];
      coords[1] = vertices[v_i[locs_hex[i][1]]];
      coords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!g_fcn_2i(metrics[i], gradients+3*i, coords, n)) return false;
    }

    switch(avgMethod) {
    case MINIMUM:
      m = metrics[0];
      for (i = 1; i < 4; ++i) {
        if (metrics[i] < m) m = metrics[i];
      }

      nm = 0;
      for (i = 0; i < 4; ++i) {
        if (metrics[i] - m <= MSQ_MIN) {
          grad[locs_hex[i][0]] += gradients[3*i+0];
          grad[locs_hex[i][1]] += gradients[3*i+1];
          grad[locs_hex[i][2]] += gradients[3*i+2];
          ++nm;
        }
      }

      for (i = 0; i < 4; ++i) {
        grad[i] /= nm;
      }
      break;

    case MAXIMUM:
      m = metrics[0];
      for (i = 1; i < 4; ++i) {
        if (metrics[i] > m) m = metrics[i];
      }

      nm = 0;
      for (i = 0; i < 4; ++i) {
        if (m - metrics[i] <= MSQ_MIN) {
          grad[locs_hex[i][0]] += gradients[3*i+0];
          grad[locs_hex[i][1]] += gradients[3*i+1];
          grad[locs_hex[i][2]] += gradients[3*i+2];
          ++nm;
        }
      }

      for (i = 0; i < 4; ++i) {
        grad[i] /= nm;
      }
      break;

    case SUM:
      m = 0;
      for (i = 0; i < 4; ++i) {
        m += metrics[i];
      }

      for (i = 0; i < 4; ++i) {
        grad[locs_hex[i][0]] += gradients[3*i+0];
        grad[locs_hex[i][1]] += gradients[3*i+1];
        grad[locs_hex[i][2]] += gradients[3*i+2];
      }
      break;

    case GEOMETRIC:
      m = 0.0;
      for (i = 0; i < 4; ++i) {
        m += log(metrics[i]);
        metrics[i] = 1.0 / metrics[i];
      }
      m = exp(m / 4.0);

      for (i = 0; i < 4; ++i) {
        grad[locs_hex[i][0]] += metrics[i]*gradients[3*i+0];
        grad[locs_hex[i][1]] += metrics[i]*gradients[3*i+1];
        grad[locs_hex[i][2]] += metrics[i]*gradients[3*i+2];
      }

      nm = m / 4.0;
      for (i = 0; i < 4; ++i) {
        grad[i] *= nm;
      }
      break;

    default:
      switch(avgMethod) {
      case LINEAR:
        t = 1.0;
        break;

      case RMS:
        t = 2.0;
        break;

      case HARMONIC:
        t = -1.0;
        break;

      case HMS:
        t = -2.0;
        break;
      default:
        err.set_msg("averaging method not available.");
        break;
      }

      m = 0;
      for (i = 0; i < 4; ++i) {
        nm = pow(metrics[i], t);
        m += nm;

        metrics[i] = t*nm/metrics[i];
      }

      nm = m / 4.0;
      m = pow(nm, 1.0 / t);

      for (i = 0; i < 4; ++i) {
        grad[locs_hex[i][0]] += metrics[i]*gradients[3*i+0];
        grad[locs_hex[i][1]] += metrics[i]*gradients[3*i+1];
        grad[locs_hex[i][2]] += metrics[i]*gradients[3*i+2];
      }

      nm = m / (4.0*nm*t);
      for (i = 0; i < 4; ++i) {
        grad[i] *= nm;
      }
      break;
    }

    // This is not very efficient, but is one way to select correct gradients
    // For gradients, info is returned only for free vertices, in the order of v[].
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < nv; ++j) {
        if (vertices + v_i[i] == v[j]) {
          g[j] = grad[i];
        }
      }
    }
    break;

  case TETRAHEDRON:
    coords[0] = vertices[v_i[0]];
    coords[1] = vertices[v_i[1]];
    coords[2] = vertices[v_i[2]];
    coords[3] = vertices[v_i[3]];
    if (!g_fcn_3e(m, grad, coords)) return false;

    // This is not very efficient, but is one way to select correct gradients.
    // For gradients, info is returned only for free vertices, in the
    // order of v[].
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < nv; ++j) {
        if (vertices + v_i[i] == v[j]) {
          g[j] = grad[i];
        }
      }
    }
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      grad[i] = 0.0;

      coords[0] = vertices[v_i[locs_hex[i][0]]];
      coords[1] = vertices[v_i[locs_hex[i][1]]];
      coords[2] = vertices[v_i[locs_hex[i][2]]];
      coords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!g_fcn_3i(metrics[i], gradients+4*i, coords)) return false;
    }

    switch(avgMethod) {
    case MINIMUM:
      m = metrics[0];
      for (i = 1; i < 8; ++i) {
        if (metrics[i] < m) m = metrics[i];
      }

      nm = 0;
      for (i = 0; i < 8; ++i) {
        if (metrics[i] - m <= MSQ_MIN) {
          grad[locs_hex[i][0]] += gradients[4*i+0];
          grad[locs_hex[i][1]] += gradients[4*i+1];
          grad[locs_hex[i][2]] += gradients[4*i+2];
          grad[locs_hex[i][3]] += gradients[4*i+3];
          ++nm;
        }
      }

      for (i = 0; i < 8; ++i) {
        grad[i] /= nm;
      }
      break;

    case MAXIMUM:
      m = metrics[0];
      for (i = 1; i < 8; ++i) {
        if (metrics[i] > m) m = metrics[i];
      }

      nm = 0;
      for (i = 0; i < 8; ++i) {
        if (m - metrics[i] <= MSQ_MIN) {
          grad[locs_hex[i][0]] += gradients[4*i+0];
          grad[locs_hex[i][1]] += gradients[4*i+1];
          grad[locs_hex[i][2]] += gradients[4*i+2];
          grad[locs_hex[i][3]] += gradients[4*i+3];
          ++nm;
        }
      }

      for (i = 0; i < 8; ++i) {
        grad[i] /= nm;
      }
      break;

    case SUM:
      m = 0;
      for (i = 0; i < 8; ++i) {
        m += metrics[i];
      }

      for (i = 0; i < 8; ++i) {
        grad[locs_hex[i][0]] += gradients[4*i+0];
        grad[locs_hex[i][1]] += gradients[4*i+1];
        grad[locs_hex[i][2]] += gradients[4*i+2];
        grad[locs_hex[i][3]] += gradients[4*i+3];
      }
      break;

    case GEOMETRIC:
      m = 0.0;
      for (i = 0; i < 8; ++i) {
        m += log(metrics[i]);
        metrics[i] = 1.0 / metrics[i];
      }
      m = exp(m / 8.0);

      for (i = 0; i < 8; ++i) {
        grad[locs_hex[i][0]] += metrics[i]*gradients[4*i+0];
        grad[locs_hex[i][1]] += metrics[i]*gradients[4*i+1];
        grad[locs_hex[i][2]] += metrics[i]*gradients[4*i+2];
        grad[locs_hex[i][3]] += metrics[i]*gradients[4*i+3];
      }

      nm = m / 8.0;
      for (i = 0; i < 8; ++i) {
        grad[i] *= nm;
      }
      break;

    default:
      switch(avgMethod) {
      case LINEAR:
        t = 1.0;
        break;

      case RMS:
        t = 2.0;
        break;

      case HARMONIC:
        t = -1.0;
        break;

      case HMS:
        t = -2.0;
        break;
      default:
        err.set_msg("averaging method not available.");
        break;
      }

      m = 0;
      for (i = 0; i < 8; ++i) {
        nm = pow(metrics[i], t);
        m += nm;

        metrics[i] = t*nm/metrics[i];
      }

      nm = m / 8.0;
      m = pow(nm, 1.0 / t);

      for (i = 0; i < 8; ++i) {
        grad[locs_hex[i][0]] += metrics[i]*gradients[4*i+0];
        grad[locs_hex[i][1]] += metrics[i]*gradients[4*i+1];
        grad[locs_hex[i][2]] += metrics[i]*gradients[4*i+2];
        grad[locs_hex[i][3]] += metrics[i]*gradients[4*i+3];
      }

      nm = m / (8.0*nm*t);
      for (i = 0; i < 8; ++i) {
        grad[i] *= nm;
      }
      break;
    }

    // This is not very efficient, but is one way to select correct gradients
    // For gradients, info is returned only for free vertices, in the order of v[].
    for (i = 0; i < 8; ++i) {
      for (j = 0; j < nv; ++j) {
        if (vertices + v_i[i] == v[j]) {
          g[j] = grad[i];
        }
      }
    }
    break;

  default:
    break;
  } // end switch over element type

//  FUNCTION_TIMER_END();
  return true;
}


#undef __FUNC__
#define __FUNC__ "MeanRatioQualityMetric::compute_element_analytical_hessian" 
bool MeanRatioQualityMetric::compute_element_analytical_hessian(PatchData &pd,
								MsqMeshEntity *e,
								MsqVertex *fv[], 
								Vector3D g[],
								Matrix3D h[],
								int nfv, 
								double &m,
								MsqError &err)
{
//  FUNCTION_TIMER_START(__FUNC__);
  EntityTopology topo = e->get_element_type();

  if (((topo == QUADRILATERAL) || (topo == HEXAHEDRON)) && 
      ((avgMethod == MINIMUM) || (avgMethod == MAXIMUM))) {
    std::cout << "Minimum and maximum not continuously differentiable." << std::endl;
    std::cout << "Element of subdifferential will be returned." << std::endl;
    std::cout << "Who knows what the Hessian is?" << std::endl;
  }

  MsqVertex *vertices = pd.get_vertex_array(err);
  const size_t *v_i = e->get_vertex_index_array();

  Vector3D coords[4];		// Vertex coordinates for the (decomposed) elements
  Vector3D n;			// Surface normal for 2D objects
  Vector3D gradients[32];	// Gradient of metric with respect to the coords
  Matrix3D hessians[80];	// Hessian of matrix with respect to the coords
  double   metrics[8];		// Metric values for the (decomposed) elements
  double   nm, t=0;

  int locs_hex[8][4] = {{0, 1, 3, 4},	// Hex element descriptions
                        {1, 2, 0, 5},
                        {2, 3, 1, 6},
                        {3, 0, 2, 7},
                        {4, 7, 5, 0},
                        {5, 4, 6, 1},
                        {6, 5, 7, 2},
                        {7, 6, 4, 3}};
  int i, j, k, l, ind;
  int r, c, loc;

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_CHKERR(err);
    n = n / n.length();		// Need unit normal
    coords[0] = vertices[v_i[0]];
    coords[1] = vertices[v_i[1]];
    coords[2] = vertices[v_i[2]];
    if (!h_fcn_2e(m, g, h, coords, n)) return false;

    // zero out fixed elements of g
    j = 0;
    for (i = 0; i < 3; ++i) {
      // if free vertex, see next
      if (vertices + v_i[i] == fv[j] )
        ++j;
      // else zero gradient entry
      else {
        g[i] = 0.;

        switch(i) {
        case 0:
          h[0] = 0.;
          h[1] = 0.;
          h[2] = 0.;
          break;

        case 1:
          h[1] = 0.;
          h[3] = 0.;
          h[4] = 0.;
          break;

        case 2:
          h[2] = 0.;
          h[4] = 0.;
          h[5] = 0.;
        }
      }
    }
    break;

  case QUADRILATERAL:
    for (i=0; i < 10; ++i) {
      h[i] = 0.;
    }
    
    pd.get_domain_normal_at_element(e, n, err); MSQ_CHKERR(err);
    for (i = 0; i < 4; ++i) {
      g[i] = 0.0;

      n = n / n.length();	// Need unit normal
      coords[0] = vertices[v_i[locs_hex[i][0]]];
      coords[1] = vertices[v_i[locs_hex[i][1]]];
      coords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!h_fcn_2i(metrics[i], gradients+3*i, hessians+6*i, coords, n)) return false;
    }

    switch(avgMethod) {
    case MINIMUM:
      err.set_msg("MINIMUM averaging method does not work.");
      return false;

    case MAXIMUM:
      err.set_msg("MAXIMUM averaging method does not work.");
      return false;

    case SUM:
      m = 0;
      for (i = 0; i < 4; ++i) {
	m += metrics[i];
      }

      l = 0;
      for (i = 0; i < 4; ++i) {
        g[locs_hex[i][0]] += gradients[3*i+0];
	g[locs_hex[i][1]] += gradients[3*i+1];
	g[locs_hex[i][2]] += gradients[3*i+2];

	for (j = 0; j < 3; ++j) {
	  for (k = j; k < 3; ++k) {
	    r = locs_hex[i][j];
	    c = locs_hex[i][k];

	    if (r <= c) {
	      loc = 4*r - (r*(r+1)/2) + c;
	      h[loc] += hessians[l];
	    } 
	    else {
	      loc = 4*c - (c*(c+1)/2) + r;
	      h[loc] += transpose(hessians[l]);
	    }
	    ++l;
	  }
	}
      }
      break;

    case GEOMETRIC:
      err.set_msg("GEOMETRIC averaging method does not work.");
      return false;

    default:
      switch(avgMethod) {
      case LINEAR:
	err.set_msg("LINEAR averaging method does not work.");
	return false;

      case RMS:
	err.set_msg("RMS averaging method does not work.");
	return false;

      case HARMONIC:
	err.set_msg("HARMONIC averaging method does not work.");
	return false;

      case HMS:
	err.set_msg("HMS averaging method does not work.");
	return false;

      default:
        err.set_msg("averaging method not available.");
        break;

      }

      m = 0;
      for (i = 0; i < 4; ++i) {
	nm = pow(metrics[i], t);
	m += nm;

	metrics[i] = t*nm/metrics[i];
      }

      nm = m / 4.0;
      m = pow(nm, 1.0 / t);

      for (i = 0; i < 4; ++i) {
        g[locs_hex[i][0]] += metrics[i]*gradients[3*i+0];
	g[locs_hex[i][1]] += metrics[i]*gradients[3*i+1];
	g[locs_hex[i][2]] += metrics[i]*gradients[3*i+2];
      }

      nm = m / (4.0*nm*t);
      for (i = 0; i < 4; ++i) {
	g[i] *= nm;
      }
      break;
    }

    // zero out fixed elements of gradient and Hessian
    ind = 0;
    for (i=0; i<4; ++i) {
      // if free vertex, see next
      if ( vertices+v_i[i] == fv[ind] )
        ++ind;
      // else zero gradient entry and hessian entries.
      else {
        g[i] = 0.;
        switch(i) {
        case 0:
          h[0]=0.;   h[1]=0.;   h[2]=0.;   h[3]=0.;
          break;
          
        case 1:
          h[1]=0.;   h[4]=0.;   h[5]=0.;   h[6]=0.;
          break;
          
        case 2:
          h[2]=0.;   h[5]=0.;   h[7]=0.;   h[8]=0.;
          break;
          
        case 3:
          h[3]=0.;   h[6]=0.;   h[8]=0.;   h[9]=0.;
          break;
        }
      }
    }

    // Makes sure we zero the hessian blocks corresponding to fixed vertices. 
    for(i=0; i<4; ++i) {
      for (j=i; j<4; ++j) {
  
        // for entry i,j in upper right part of 4*4 matrix, index in a 1D array is
        ind = 4*i- (i*(i+1)/2) +j; // 4*i - \sum_{n=0}^i n +j 

        bool nul_hessian = true; 
        for (k=0; k<nfv; ++k) 
          for (l=0; l<nfv; ++l) 
            if ( (vertices + v_i[i] == fv[k]) && (vertices + v_i[j] == fv[l]) )
              nul_hessian = false;

        if (nul_hessian == true)
          h[ind] = 0.;
      }
    }
    break;

  case TETRAHEDRON:
    coords[0] = vertices[v_i[0]];
    coords[1] = vertices[v_i[1]];
    coords[2] = vertices[v_i[2]];
    coords[3] = vertices[v_i[3]];
    if (!h_fcn_3e(m, g, h, coords)) return false;

    // zero out fixed elements of g
    j = 0;
    for (i = 0; i < 4; ++i) {
      // if free vertex, see next
      if (vertices + v_i[i] == fv[j] )
        ++j;
      // else zero gradient entry
      else {
        g[i] = 0.;

        switch(i) {
        case 0:
          h[0] = 0.;
          h[1] = 0.;
          h[2] = 0.;
          h[3] = 0.;
          break;

        case 1:
          h[1] = 0.;
          h[4] = 0.;
          h[5] = 0.;
          h[6] = 0.;
          break;

        case 2:
          h[2] = 0.;
          h[5] = 0.;
          h[7] = 0.;
          h[8] = 0.;
          break;

        case 3:
          h[3] = 0.;
          h[6] = 0.;
          h[8] = 0.;
          h[9] = 0.;
          break;
        }
      }
    }
    break;

  case HEXAHEDRON:
    for (i=0; i<36; ++i)
      h[i] = 0.;

    for (i = 0; i < 8; ++i) {
      g[i] = 0.0;

      coords[0] = vertices[v_i[locs_hex[i][0]]];
      coords[1] = vertices[v_i[locs_hex[i][1]]];
      coords[2] = vertices[v_i[locs_hex[i][2]]];
      coords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!h_fcn_3i(metrics[i], gradients+4*i, hessians+10*i, coords)) return false;
    }

    switch(avgMethod) {
    case MINIMUM:
      err.set_msg("MINIMUM averaging method does not work.");
      return false;
//       m = metrics[0];
//       for (i = 1; i < 8; ++i) {
// 	if (metrics[i] < m) m = metrics[i];
//       }

//       nm = 0;
//       for (i = 0; i < 8; ++i) {
//         if (metrics[i] <= m + MSQ_MIN) {
// 	  g[locs_hex[i][0]] += gradients[4*i+0];
// 	  g[locs_hex[i][1]] += gradients[4*i+1];
// 	  g[locs_hex[i][2]] += gradients[4*i+2];
// 	  g[locs_hex[i][3]] += gradients[4*i+3];
// 	  ++nm;
//         }
//       }

//       for (i = 0; i < 8; ++i) {
// 	g[i] /= nm;
//       }
//       break;

    case MAXIMUM:
      err.set_msg("MAXIMUM averaging method does not work.");
      return false;
//       m = metrics[0];
//       for (i = 1; i < 8; ++i) {
// 	if (metrics[i] > m) m = metrics[i];
//       }

//       nm = 0;
//       for (i = 0; i < 8; ++i) {
//         if (metrics[i] >= m - MSQ_MIN) {
// 	  g[locs_hex[i][0]] += gradients[4*i+0];
// 	  g[locs_hex[i][1]] += gradients[4*i+1];
// 	  g[locs_hex[i][2]] += gradients[4*i+2];
// 	  g[locs_hex[i][3]] += gradients[4*i+3];
// 	  ++nm;
//         }
//       }

//       for (i = 0; i < 8; ++i) {
// 	g[i] /= nm;
//       }
//       break;

    case SUM:
      m = 0;
      for (i = 0; i < 8; ++i) {
	m += metrics[i];
      }

      l = 0;
      for (i = 0; i < 8; ++i) {
        g[locs_hex[i][0]] += gradients[4*i+0];
	g[locs_hex[i][1]] += gradients[4*i+1];
	g[locs_hex[i][2]] += gradients[4*i+2];
	g[locs_hex[i][3]] += gradients[4*i+3];

	for (j = 0; j < 4; ++j) {
	  for (k = j; k < 4; ++k) {
	    r = locs_hex[i][j];
	    c = locs_hex[i][k];

	    if (r <= c) {
	      loc = 8*r - (r*(r+1)/2) + c;
	      h[loc] += hessians[l];
	    } 
	    else {
	      loc = 8*c - (c*(c+1)/2) + r;
	      h[loc] += transpose(hessians[l]);
	    }
	    ++l;
	  }
	}
      }
      break;

    case GEOMETRIC:
      err.set_msg("GEOMETRIC averaging method does not work.");
      return false;
//       m = 0.0;
//       for (i = 0; i < 8; ++i) {
// 	m += log(metrics[i]);
// 	metrics[i] = 1.0 / metrics[i];
//       }
//       m = exp(m / 8.0);

//       for (i = 0; i < 8; ++i) {
//         g[locs_hex[i][0]] += metrics[i]*gradients[4*i+0];
// 	g[locs_hex[i][1]] += metrics[i]*gradients[4*i+1];
// 	g[locs_hex[i][2]] += metrics[i]*gradients[4*i+2];
// 	g[locs_hex[i][3]] += metrics[i]*gradients[4*i+3];
//       }

//       nm = m / 8.0;
//       for (i = 0; i < 8; ++i) {
// 	g[i] *= nm;
//       }
//       break;

    default:
      switch(avgMethod) {
      case LINEAR:
	err.set_msg("LINEAR averaging method does not work.");
	return false;
// 	t = 1.0;
// 	break;

      case RMS:
	err.set_msg("RMS averaging method does not work.");
	return false;
// 	t = 2.0;
// 	break;

      case HARMONIC:
	err.set_msg("HARMONIC averaging method does not work.");
	return false;
// 	t = -1.0;
// 	break;

      case HMS:
	err.set_msg("HMS averaging method does not work.");
	return false;
// 	t = -2.0;
// 	break;
      default:
        err.set_msg("averaging method not available.");
        break;

      }

      m = 0;
      for (i = 0; i < 8; ++i) {
	nm = pow(metrics[i], t);
	m += nm;

	metrics[i] = t*nm/metrics[i];
      }

      nm = m / 8.0;
      m = pow(nm, 1.0 / t);

      for (i = 0; i < 8; ++i) {
        g[locs_hex[i][0]] += metrics[i]*gradients[4*i+0];
	g[locs_hex[i][1]] += metrics[i]*gradients[4*i+1];
	g[locs_hex[i][2]] += metrics[i]*gradients[4*i+2];
	g[locs_hex[i][3]] += metrics[i]*gradients[4*i+3];
      }

      nm = m / (8.0*nm*t);
      for (i = 0; i < 8; ++i) {
	g[i] *= nm;
      }
      break;
    }

    // zero out fixed elements of gradient and Hessian
    ind = 0;
    for (i=0; i<8; ++i) {
      // if free vertex, see next
      if ( vertices+v_i[i] == fv[ind] )
        ++ind;
      // else zero gradient entry and hessian entries.
      else {
        g[i] = 0.;
        switch(i) {
        case 0:
          h[0]=0.;   h[1]=0.;   h[2]=0.;   h[3]=0.;
          h[4]=0.;   h[5]=0.;   h[6]=0.;   h[7]=0.;
          break;
          
        case 1:
          h[1]=0.;   h[8]=0.;   h[9]=0.;   h[10]=0.;
          h[11]=0.;  h[12]=0.;  h[13]=0.;  h[14]=0.;
          break;
          
        case 2:
          h[2]=0.;   h[9]=0.;   h[15]=0.;  h[16]=0.;
          h[17]=0.;  h[18]=0.;  h[19]=0.;  h[20]=0.;
          break;
          
        case 3:
          h[3]=0.;   h[10]=0.;  h[16]=0.;  h[21]=0.;
          h[22]=0.;  h[23]=0.;  h[24]=0.;  h[25]=0.;
          break;
          
        case 4:
          h[4]=0.;   h[11]=0.;  h[17]=0.;  h[22]=0.;
          h[26]=0.;  h[27]=0.;  h[28]=0.;  h[29]=0.;
          break;
          
        case 5:
          h[5]=0.;   h[12]=0.;  h[18]=0.;  h[23]=0.;
          h[27]=0.;  h[30]=0.;  h[31]=0.;  h[32]=0.;
          break;
          
        case 6:
          h[6]=0.;   h[13]=0.;  h[19]=0.;  h[24]=0.;
          h[28]=0.;  h[31]=0.;  h[33]=0.;  h[34]=0.;
          break;
          
        case 7:
          h[7]=0.;   h[14]=0.;  h[20]=0.;  h[25]=0.;
          h[29]=0.;  h[32]=0.;  h[34]=0.;  h[35]=0.;
          break;
        }
      }
    }

    // Makes sure we zero the hessian blocks corresponding to fixed vertices. 
    for(i=0; i<8; ++i) {
      for (j=i; j<8; ++j) {
  
        // for entry i,j in upper right part of 8*8 matrix, index in a 1D array is
        ind = 8*i- (i*(i+1)/2) +j; // 8*i - \sum_{n=0}^i n +j 

        bool nul_hessian = true; 
        for (k=0; k<nfv; ++k) 
          for (l=0; l<nfv; ++l) 
            if ( (vertices + v_i[i] == fv[k]) && (vertices + v_i[j] == fv[l]) )
              nul_hessian = false;

        if (nul_hessian == true)
          h[ind] = 0.;
      }
    }

    break;

  default:
    break;
  } // end switch over element type

//  FUNCTION_TIMER_END();
  return true;
}

