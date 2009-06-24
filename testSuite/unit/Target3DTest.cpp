/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Target2DTest.cpp
 *  \brief Unit tests for 2D target metrics
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TargetMetric3D.hpp"
#include "UnitUtil.hpp"
#include "MsqError.hpp"

using namespace Mesquite;

// Test functions implemented in class TargetMetric3D
class TargetMetric3DTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TargetMetric3DTest );
  CPPUNIT_TEST (test_numerical_gradient);
  CPPUNIT_TEST (test_numerical_hessian);
  CPPUNIT_TEST_SUITE_END(); 
  public:
  void test_numerical_gradient();
  void test_numerical_hessian();
};

// Test concrete implementations of TargetMetric3D
template <class Metric>
class Target3DTest : public CppUnit::TestFixture
{
private:
  Metric test_metric;
  bool shapeInvariant, sizeInvariant, orientInvariant, Barrier;
  double idealVal;
public:
  Target3DTest( bool shape_invariant, 
                bool size_invariant, 
                bool orient_invariant, 
                bool barrier, 
                double ideal_element_val )
    : shapeInvariant(shape_invariant),
      sizeInvariant(size_invariant), 
      orientInvariant(orient_invariant), 
      Barrier(barrier),
      idealVal(ideal_element_val)
    {}
  
  void test_ideal_eval();
  void test_ideal_gradient();
  void test_inverted();
  void test_shape();
  void test_scale();
  void test_orient();
  
  void compare_anaytic_and_numeric_grads();
  void compare_anaytic_and_numeric_hess();
  void compare_eval_and_eval_with_grad();
  void compare_eval_with_grad_and_eval_with_hess();
  
  void test_non_ideal( bool sensitive,
                       MsqMatrix<3,3> A,
                       MsqMatrix<3,3> B );
};
  
template <class Metric> void Target3DTest<Metric>::test_ideal_eval()
{
  MsqPrintError err(msq_stdio::cerr);
  const double Avals[] = { 2, 1, 1, 1, 2, 1, 1, 1, 2 };
  const double Bvals[] = { 1.5, -0.7, -0.8, 0.8, -1.3, -0.7, 0.6, -0.9, -2.0 };
  const MsqMatrix<3,3> I(1.0), A(Avals), B(Bvals);
  double val, eps = 1e-6;;
  bool valid;
  
  valid = test_metric.evaluate( I, I, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
  
  valid = test_metric.evaluate( A, A, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
  
  valid = test_metric.evaluate( B, B, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
}
  
template <class Metric> void Target3DTest<Metric>::test_ideal_gradient()
{
  MsqPrintError err(msq_stdio::cerr);
  const double Avals[] = { 2, 1, 1, 1, 2, 1, 1, 1, 2 };
  const double Bvals[] = { 1.5, -0.7, -0.8, 0.8, -1.3, -0.7, 0.6, -0.9, -2.0 };
  const MsqMatrix<3,3> I(1.0), A(Avals), B(Bvals);
  MsqMatrix<3,3> grad;
  double val, eps = 1e-6;
  bool valid;
  
  valid = test_metric.evaluate_with_grad( I, I, val, grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(0.0)), grad, eps );
  
  valid = test_metric.evaluate_with_grad( A, A, val, grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(0.0)), grad, eps );
  
  valid = test_metric.evaluate_with_grad( B, B, val, grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(0.0)), grad, eps );
}

template <class Metric> void Target3DTest<Metric>::test_inverted() 
{
  MsqPrintError err(msq_stdio::cerr);
  const double A_vals[] = { 1, 0, 0, 
                            0, 1, 0, 
                            0, 0, -1 };
  MsqMatrix<3,3> A( A_vals ), W( 1.0 ), grad, hess[6];
  double val;
  bool valid;
  
  if (Barrier) {
    valid = test_metric.evaluate( A, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!valid);
    
    valid = test_metric.evaluate_with_grad( A, W, val, grad, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!valid);
    
    valid = test_metric.evaluate_with_hess( A, W, val, grad, hess, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!valid);
  }
  else {
    valid = test_metric.evaluate( A, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT( val > idealVal );
    
    valid = test_metric.evaluate_with_grad( A, W, val, grad, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT( val > idealVal );
  }
}
  
template <class Metric> void Target3DTest<Metric>::test_non_ideal( 
                                                     bool sensitive,
                                                     MsqMatrix<3,3> A,
                                                     MsqMatrix<3,3> W )
{
  MsqPrintError err(msq_stdio::cerr);
  MsqMatrix<3,3> grad;
  double val, eps = 1e-6;
  bool valid;
  if (!sensitive) {
    valid = test_metric.evaluate( A, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );

    valid = test_metric.evaluate_with_grad( A, W, val, grad, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
    ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(0.0)), grad, eps );
  }
  else {
    valid = test_metric.evaluate( A, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT( val > idealVal );
  }
}

template <class Metric> void Target3DTest<Metric>::test_shape()
{
  const double r3 = sqrt(3.0);
  const double A_vals[] = { 2/r3, 1/r3, 0, 
                            1/r3, 2/r3, 0, 
                            0,    0,    1 };
  MsqMatrix<3,3> A( A_vals ), W( 1.0 );
  test_non_ideal( !shapeInvariant, A, W );
}
  
template <class Metric> void Target3DTest<Metric>::test_scale() 
{
  MsqMatrix<3,3> A( 2.0 ), W( 1.0 );
  test_non_ideal( !sizeInvariant, A, W );
}
  
template <class Metric> void Target3DTest<Metric>::test_orient() 
{
  const double A_vals[] = { 0, -1,  0, 
                            1,  0,  0, 
                            0,  0,  1 };
  MsqMatrix<3,3> A( A_vals ), W( 1.0 );
  test_non_ideal( !orientInvariant, A, W );
}

template <class Metric> 
void Target3DTest<Metric>::compare_eval_and_eval_with_grad()
{
  Metric metric;
  
  const double Avals[] = { 2, 1, 1, 1, 2, 1, 1, 1, 2 };
  const double Bvals[] = { 1.5, -0.7, -0.8, 0.8, -1.3, -0.7, 0.6, -0.9, -2.0 };
  const MsqMatrix<3,3> I( 1.0 );
  const MsqMatrix<3,3> A( Avals );
  const MsqMatrix<3,3> B( Bvals );
  
  MsqError err;
  MsqMatrix<3,3> g;
  bool valid;
  double gv, v;
  
  valid = metric.evaluate_with_grad( I, A, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate( I, A, v, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v, gv, 1e-6 );
  
  valid = metric.evaluate_with_grad( A, B, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate( A, B, v, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v, gv, 1e-6 );
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  const double Cvals[] = { 0.5, 0.0,  0.1,
                           0.5, 1.0,  0.1,
                           0.0, 0.0, -1.5 };
  const MsqMatrix<3,3> C(Cvals);
  valid = metric.evaluate_with_grad( C, I, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate( C, I, v, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v, gv, 1e-6 );
}

template <class Metric> 
void Target3DTest<Metric>::compare_eval_with_grad_and_eval_with_hess()
{
  Metric metric;
  
  const double Avals[] = { 2, 1, 1, 1, 2, 1, 1, 1, 2 };
  const double Bvals[] = { 1.5, -0.7, -0.8, 0.8, -1.3, -0.7, 0.6, -0.9, -2.0 };
  const MsqMatrix<3,3> I( 1.0 );
  const MsqMatrix<3,3> A( Avals );
  const MsqMatrix<3,3> B( Bvals );
  
  MsqError err;
  MsqMatrix<3,3> g, h, hess[6];
  bool valid;
  double gv, hv;
  
  valid = metric.evaluate_with_grad( I, A, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, A, hv, h, hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( gv, hv, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, h, 1e-5 );
  
  valid = metric.evaluate_with_grad( A, B, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( A, B, hv, h, hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( gv, hv, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, h, 1e-5 );
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  const double Cvals[] = { 0.5, 0.0,  0.1,
                           0.5, 1.0,  0.1,
                           0.0, 0.0, -1.5 };
  const MsqMatrix<3,3> C(Cvals);
  valid = metric.evaluate_with_grad( C, I, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( C, I, hv, h, hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( gv, hv, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, h, 1e-5 );
}

template <class Metric> 
void Target3DTest<Metric>::compare_anaytic_and_numeric_grads()
{
  const double EPS_VAL = 1e-6;
  const double EPS_GRAD = 5e-3;
  Metric metric;
  
  const double Avals[] = { 2, 1, 1, 1, 2, 1, 1, 1, 2 };
  const double Bvals[] = { 1.5, -0.7, -0.8, 0.8, -1.3, -0.7, 0.6, -0.9, -2.0 };
  const MsqMatrix<3,3> I( 1.0 );
  const MsqMatrix<3,3> A( Avals );
  const MsqMatrix<3,3> B( Bvals );
  
  MsqError err;
  MsqMatrix<3,3> num, ana;
  bool valid;
  double nval, aval;
  
  MsqMatrix<3,3> D(I);
  D(0,0) += 1e-5;
  valid = metric.TargetMetric3D::evaluate_with_grad( D, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( D, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
  
  valid = metric.TargetMetric3D::evaluate_with_grad( I, A, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( I, A, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
  
  valid = metric.TargetMetric3D::evaluate_with_grad( A, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
  
  valid = metric.TargetMetric3D::evaluate_with_grad( I, B, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( I, B, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
  
  valid = metric.TargetMetric3D::evaluate_with_grad( B, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( B, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
   
  valid = metric.TargetMetric3D::evaluate_with_grad( A, B, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, B, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
  
  valid = metric.TargetMetric3D::evaluate_with_grad( A, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  const double Cvals[] = { 0.5, 0.0,  0.1,
                           0.5, 1.0,  0.1,
                           0.0, 0.0, -1.5 };
  const MsqMatrix<3,3> C(Cvals);
  valid = metric.TargetMetric3D::evaluate_with_grad( C, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( C, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, EPS_GRAD );
}


template <class Metric> 
void Target3DTest<Metric>::compare_anaytic_and_numeric_hess()
{
  const double EPS_VAL = 1e-6;
  const double EPS_GRAD = 5e-3;
  const double EPS_HESS = EPS_GRAD;
  Metric metric;
  
  const double Avals[] = { 2, 1, 1, 1, 2, 1, 1, 1, 2 };
  const double Bvals[] = { 1.5, -0.7, -0.8, 0.8, -1.3, -0.7, 0.6, -0.9, -2.0 };
  const MsqMatrix<3,3> I( 1.0 );
  const MsqMatrix<3,3> A( Avals );
  const MsqMatrix<3,3> B( Bvals );
  
  MsqError err;
  MsqMatrix<3,3> dmdA_num, dmdA_ana, d2mdA2_num[6], d2mdA2_ana[6];
  bool valid;
  double val_num, val_ana;
  
  valid = metric.TargetMetric3D::evaluate_with_hess( I, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
  
  valid = metric.TargetMetric3D::evaluate_with_hess( I, A, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, A, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
  
  valid = metric.TargetMetric3D::evaluate_with_hess( A, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( A, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
 
  valid = metric.TargetMetric3D::evaluate_with_hess( B, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( B, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
  
  valid = metric.TargetMetric3D::evaluate_with_hess( I, B, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, B, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
  
  valid = metric.TargetMetric3D::evaluate_with_hess( A, B, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( A, B, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
  
  valid = metric.TargetMetric3D::evaluate_with_hess( B, A, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( B, A, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  const double Cvals[] = { 0.5, 0.0,  0.1,
                           0.5, 1.0,  0.1,
                           0.0, 0.0, -1.5 };
  const MsqMatrix<3,3> C(Cvals);
  valid = metric.TargetMetric3D::evaluate_with_hess( C, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( C, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, EPS_GRAD );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], EPS_HESS );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], EPS_HESS );
}

// implement metric that is the sum of the elements of 2A-W,
// such that the derivative of the result with resepct to
// each element of A is 2.
class GradTestMetric3D : public TargetMetric3D
{
  public:
    msq_std::string get_name() const { return "GradTest"; }
  
    bool evaluate( const MsqMatrix<3,3>& A,
                   const MsqMatrix<3,3>& W,
                   double& result,
                   MsqError&  )
    {
      result = 0;
      for (int r = 0; r < 3; ++r) 
        for (int c = 0; c < 3; ++c)
          result += 2 * A(r,c) - W(r,c);
      return true;
    }
};

// implement metric: |2A - A^T - W|^2
// such that the Hessian is the constant:
//  _                         _
// | 2  0  0  0  0  0  0  0  0 |
// | 0 10  0 -8  0  0  0  0  0 |
// | 0  0 10  0  0  0 -8  0  0 |
// | 0 -8  0 10  0  0  0  0  0 |
// | 0  0  0  0  2  0  0  0  0 |
// | 0  0  0  0  0 10  0 -8  0 |
// | 0  0 -8  0  0  0 10  0  0 |
// | 0  0  0  0  0 -8  0 10  0 |
// |_0  0  0  0  0  0  0  0  2_|

class HessTestMetric3D : public TargetMetric3D
{
  public:
    msq_std::string get_name() const { return "HessTest"; }

    bool evaluate( const MsqMatrix<3,3>& A,
                   const MsqMatrix<3,3>& W,
                   double& result,
                   MsqError&  )
    {
      result = sqr_Frobenius(2*A - transpose(A) - W);
      return true;
    }
    bool evaluate_with_grad( const MsqMatrix<3,3>& A,
                             const MsqMatrix<3,3>& W,
                             double& result,
                             MsqMatrix<3,3>& wrt_A,
                             MsqError&  )
    {
      result = sqr_Frobenius(2*A - transpose(A) - W);
      wrt_A = 10*A - 8*transpose(A);
      return true;
    }
};

/** Simple target metric for testing second partial derivatives.  
 *  \f$\mu(A,W) = |A|\f$
 *  \f$\frac{\partial\mu}{\partial A} = \frac{1}{|A|} A \f$
 *  \f$\frac{\partial^{2}\mu}{\partial a_{i,i}^2} = \frac{1}{|A|} - \frac{a_{i,i}^2}{|A|^3}\f$
 *  \f$\frac{\partial^{2}\mu}{\partial a_{i,j} \partial a_{k,l} (i \ne k or j \ne l)} = -\frac{a_{i,j} a_{k,l}}{|A|^3}\f$
 */
class HessTestMetric3D_2 : public TargetMetric3D
{
  public:
    msq_std::string get_name() const { return "HessTest2"; }
  
    bool evaluate( const MsqMatrix<3,3>& A, const MsqMatrix<3,3>&, double& result, MsqError& err )
      { result = Frobenius(A); return true; }
    
    bool evaluate_with_grad( const MsqMatrix<3,3>& A, 
                             const MsqMatrix<3,3>&,
                             double& result,
                             MsqMatrix<3,3>& d,
                             MsqError& err )
    {
      result = Frobenius(A);
      d = A / result;
      return true;
    }
    
    bool evaluate_with_hess( const MsqMatrix<3,3>& A, 
                             const MsqMatrix<3,3>&,
                             double& result,
                             MsqMatrix<3,3>& d,
                             MsqMatrix<3,3> d2[6],
                             MsqError& err )
    {
      result = Frobenius(A);
      d = A / result;
      int h = 0;
      for (int r = 0; r < 3; ++r) {
        int i = h;
        for (int c = r; c < 3; ++c)
          d2[h++] = transpose(A.row(r)) * A.row(c) / -(result*result*result);
        d2[i] += MsqMatrix<3,3>(1.0/result);
      }    
      return true;
    }
};

void TargetMetric3DTest::test_numerical_gradient()
{
  GradTestMetric3D metric;
  HessTestMetric3D_2 metric2;
  const double Avals[] = { 1, 2, 3, 4, 1, 4, 3, 2, 1 };
  const double Bvals[] = { 0.1, 0.15, 0.05, 0.2, -0.1, -0.15, -0.05, -0.2, 2 };
  const MsqMatrix<3,3> I( 1.0 );
  const MsqMatrix<3,3> A( Avals );
  const MsqMatrix<3,3> B( Bvals );
  
  MsqError err;
  MsqMatrix<3,3> d;
  bool valid;
  double val, gval;
  
  valid = metric.evaluate( I, A, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( I, A, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,2), 1e-6 );
  
  valid = metric.evaluate( A, I, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, I, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,2), 1e-6 );
  
  valid = metric.evaluate( I, B, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( I, B, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,2), 1e-6 );
  
  valid = metric.evaluate( B, I, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( B, I, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,2), 1e-6 );
   
  valid = metric.evaluate( A, B, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, B, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,2), 1e-6 );
  
  valid = metric.evaluate( B, A, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( B, A, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(2,2), 1e-6 );
  
  MsqMatrix<3,3> da;
  valid = metric.evaluate_with_grad( A, I, val, da, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.TargetMetric3D::evaluate_with_grad( A, I, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  ASSERT_MATRICES_EQUAL( da, d, 1e-6 );
  
  valid = metric.evaluate_with_grad( B, I, val, da, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.TargetMetric3D::evaluate_with_grad( B, I, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  ASSERT_MATRICES_EQUAL( da, d, 1e-6 );
}


void TargetMetric3DTest::test_numerical_hessian()
{
  HessTestMetric3D metric;
  HessTestMetric3D_2 metric2;
  const double Avals[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  const double Bvals[] = { 0.1, 0.15, 0.05, 0.2, -0.1, -0.15, -0.05, -0.2, 2 };
  const MsqMatrix<3,3> I( 1.0 );
  const MsqMatrix<3,3> A( Avals );
  const MsqMatrix<3,3> B( Bvals );
  
  MsqError err;
  MsqMatrix<3,3> g, gh;
  MsqMatrix<3,3> h[6];
  bool valid;
  double val, hval;
  
  const double h_00[] = { 2, 0, 0, 
                          0,10, 0,
                          0, 0,10};
  const double h_01[] = { 0, 0, 0, 
                         -8, 0, 0,
                          0, 0, 0};
  const double h_02[] = { 0, 0, 0, 
                          0, 0, 0,
                         -8, 0, 0};
  const double h_11[] = {10, 0, 0, 
                          0, 2, 0,
                          0, 0,10};
  const double h_12[] = { 0, 0, 0, 
                          0, 0, 0,
                          0,-8, 0};
  const double h_22[] = {10, 0, 0, 
                          0,10, 0,
                          0, 0, 2};
  MsqMatrix<3,3> h00(h_00), h01(h_01), h02(h_02), h11(h_11), h12(h_12), h22(h_22);
  
  valid = metric.evaluate_with_grad( I, A, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, A, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( h00, h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( h01, h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( h02, h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( h11, h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( h12, h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( h22, h[5], 1e-6 );
  
  valid = metric.evaluate_with_grad( A, I, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( A, I, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( h00, h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( h01, h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( h02, h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( h11, h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( h12, h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( h22, h[5], 1e-6 );
  
  valid = metric.evaluate_with_grad( I, B, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, B, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( h00, h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( h01, h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( h02, h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( h11, h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( h12, h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( h22, h[5], 1e-6 );
  
  valid = metric.evaluate_with_grad( B, I, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( B, I, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( h00, h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( h01, h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( h02, h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( h11, h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( h12, h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( h22, h[5], 1e-6 );
  
  valid = metric.evaluate_with_grad( A, B, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( A, B, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( h00, h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( h01, h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( h02, h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( h11, h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( h12, h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( h22, h[5], 1e-6 );
  
  valid = metric.evaluate_with_grad( B, A, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( B, A, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( h00, h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( h01, h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( h02, h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( h11, h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( h12, h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( h22, h[5], 1e-6 );
  
  MsqMatrix<3,3> ah[6];
  valid = metric2.evaluate_with_hess( A, I, val, g, ah, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric2.TargetMetric3D::evaluate_with_hess( A, I, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[0], h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[1], h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[2], h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[3], h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[4], h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[5], h[5], 1e-6 );

  valid = metric2.evaluate_with_hess( B, I, val, g, ah, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric2.TargetMetric3D::evaluate_with_hess( B, I, hval, gh, h, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, hval, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, gh, 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[0], h[0], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[1], h[1], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[2], h[2], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[3], h[3], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[4], h[4], 1e-6 );
  ASSERT_MATRICES_EQUAL( ah[5], h[5], 1e-6 );
}


#include "TSquared3D.hpp"
#include "Target3DShapeSizeOrient.hpp"
#include "Target3DShapeSize.hpp"
#include "Target3DShapeBarrier.hpp"
#include "Target3DShape.hpp"
#include "Target3DShapeSizeBarrierAlt1.hpp"
#include "Target3DShapeSizeBarrierAlt2.hpp"
#include "InverseMeanRatio3D.hpp"
#include "Target3DSize.hpp"
#include "Target3DSizeBarrier.hpp"
#include "Target3DShapeSizeOrientBarrier.hpp"
#include "Target3DShapeSizeOrientBarrierAlt2.hpp"
#include "Target3DShapeOrientAlt1.hpp"
#include "Target3DShapeOrientBarrierAlt1.hpp"
#include "Target3DUntangle.hpp"
#include "Target3DUntangleAlt1.hpp"

#define REGISTER_TARGET3D_TEST( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class Test_ ## METRIC : public Target3DTest<METRIC> { public: \
  Test_ ## METRIC () : Target3DTest<METRIC>( (SHAPE_INVAR), (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( Test_ ## METRIC ); \
  CPPUNIT_TEST (test_ideal_eval); \
  CPPUNIT_TEST (test_ideal_gradient); \
  CPPUNIT_TEST (test_inverted); \
  CPPUNIT_TEST (test_shape); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _FileRegister ("Target3DTest"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _BaseRegister ( "Test_" #METRIC )

#define REGISTER_TARGET3D_TEST_WITH_GRAD( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class Test_ ## METRIC : public Target3DTest<METRIC> { public: \
  Test_ ## METRIC () : Target3DTest<METRIC>( (SHAPE_INVAR), (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( Test_ ## METRIC ); \
  CPPUNIT_TEST (test_ideal_eval); \
  CPPUNIT_TEST (test_ideal_gradient); \
  CPPUNIT_TEST (test_inverted); \
  CPPUNIT_TEST (test_shape); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST (compare_eval_and_eval_with_grad); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_grads); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _FileRegister ("Target3DTest"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _BaseRegister ( "Test_" #METRIC )

#define REGISTER_TARGET3D_TEST_WITH_HESS( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class Test_ ## METRIC : public Target3DTest<METRIC> { public: \
  Test_ ## METRIC () : Target3DTest<METRIC>( (SHAPE_INVAR), (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( Test_ ## METRIC ); \
  CPPUNIT_TEST (test_ideal_eval); \
  CPPUNIT_TEST (test_ideal_gradient); \
  CPPUNIT_TEST (test_inverted); \
  CPPUNIT_TEST (test_shape); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST (compare_eval_and_eval_with_grad); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_grads); \
  CPPUNIT_TEST (compare_eval_with_grad_and_eval_with_hess); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_hess); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _FileRegister ("Target3DTest"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _BaseRegister ( "Test_" #METRIC )

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric3DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric3DTest, "Target3DTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric3DTest, "TargetMetric3DTest" );
//                                Metric                             !shape !size !orient barrer ideal
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeSizeOrient,           false, false, false, false, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeSize,                 false, false,  true,  true, 1.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShape,                     false,  true,  true, false, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeBarrier,              false,  true,  true,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeOrientAlt1,           false,  true, false, false, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeOrientBarrierAlt1,    false,  true, false,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeSizeBarrierAlt1,      false, false,  true,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeSizeBarrierAlt2,      false, false,  true,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( InverseMeanRatio3D,                false,  true,  true,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DSize,                       true, false,  true, false, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DSizeBarrier,                true, false,  true,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeSizeOrientBarrier,    false, false, false,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DShapeSizeOrientBarrierAlt2,false, false, false,  true, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DUntangle,                   true,  true,  true, false, 0.0 );
REGISTER_TARGET3D_TEST_WITH_HESS( Target3DUntangleAlt1,               true,  true,  true, false, 0.0 );


class Test_TSquared3D : public Target3DTest<TSquared3D> {
  public: 
    Test_TSquared3D() : Target3DTest<TSquared3D>(false,false,false,false,0.0) {}
    CPPUNIT_TEST_SUITE( Test_TSquared3D );
    CPPUNIT_TEST( compare_eval_and_eval_with_grad ); 
    CPPUNIT_TEST( compare_anaytic_and_numeric_grads );
    CPPUNIT_TEST( compare_eval_with_grad_and_eval_with_hess );
    CPPUNIT_TEST( compare_anaytic_and_numeric_hess );
    CPPUNIT_TEST_SUITE_END();
};
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( Test_TSquared3D, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( Test_TSquared3D, "Target3DTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( Test_TSquared3D, "Test_TSquared3D" );
