/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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

    (2006) kraftche@cae.wisc.edu    
    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Target2DTest.cpp
 *  \brief Unit tests for 2D target metrics
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TargetMetric2D.hpp"
#include "UnitUtil.hpp"
#include "MsqError.hpp"

using namespace Mesquite;

class TargetMetric2DTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TargetMetric2DTest );
  CPPUNIT_TEST (test_numerical_gradient);
  CPPUNIT_TEST (test_numerical_hessian);
  CPPUNIT_TEST_SUITE_END(); 
  public:
  void test_numerical_gradient();
  void test_numerical_hessian();
};

template <class Metric>
class Target2DTest : public CppUnit::TestFixture
{
private:
  Metric test_metric;
  bool shapeInvariant, sizeInvariant, orientInvariant, Barrier;
  double idealVal;
public:
  Target2DTest( bool shape_invariant, 
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
                       MsqMatrix<2,2> A,
                       MsqMatrix<2,2> B );
};
  
template <class Metric> void Target2DTest<Metric>::test_ideal_eval()
{
  MsqPrintError err(std::cerr);
  const double Avals[] = { 2, 1, 1, 2 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I(1.0), A(Avals), B(Bvals);
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
template <class Metric> void Target2DTest<Metric>::test_ideal_gradient()
{
  MsqPrintError err(std::cerr);
  const double Avals[] = { 2, 1, 1, 2 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I(1.0), A(Avals), B(Bvals);
  MsqMatrix<2,2> grad;
  double val, eps = 5e-3;
  bool valid;
  
  valid = test_metric.evaluate_with_grad( I, I, val, grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( (MsqMatrix<2,2>(0.0)), grad, eps );
  
  valid = test_metric.evaluate_with_grad( A, A, val, grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( (MsqMatrix<2,2>(0.0)), grad, eps );
  
  valid = test_metric.evaluate_with_grad( B, B, val, grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( (MsqMatrix<2,2>(0.0)), grad, eps );
}

template <class Metric> void Target2DTest<Metric>::test_inverted() 
{
  MsqPrintError err(std::cerr);
  const double A_vals[] = { 1,  0, 
                            0, -1 };
  MsqMatrix<2,2> A( A_vals ), W( 1.0 ), grad, hess[3];
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
  
template <class Metric> void Target2DTest<Metric>::test_non_ideal( 
                                                     bool sensitive,
                                                     MsqMatrix<2,2> A,
                                                     MsqMatrix<2,2> W )
{
  MsqPrintError err(std::cerr);
  MsqMatrix<2,2> grad;
  double val, eps = 1e-5;
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
    ASSERT_MATRICES_EQUAL( (MsqMatrix<2,2>(0.0)), grad, eps );
  }
  else {
    valid = test_metric.evaluate( A, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT( val > idealVal );
  }
}

template <class Metric> void Target2DTest<Metric>::test_shape()
{
  const double r3 = sqrt(3.0);
  const double A_vals[] = { 2/r3, 1/r3, 
                            1/r3, 2/r3 };
  MsqMatrix<2,2> A( A_vals ), W( 1.0 );
  test_non_ideal( !shapeInvariant, A, W );
}
  
template <class Metric> void Target2DTest<Metric>::test_scale() 
{
  MsqMatrix<2,2> A( 2.0 ), W( 1.0 );
  test_non_ideal( !sizeInvariant, A, W );
}
  
template <class Metric> void Target2DTest<Metric>::test_orient() 
{
  const double A_vals[] = { 0, -1, 
                            1,  0 };
  MsqMatrix<2,2> A( A_vals ), W( 1.0 );
  test_non_ideal( !orientInvariant, A, W );
}


template <class Metric> 
void Target2DTest<Metric>::compare_eval_and_eval_with_grad()
{
  Metric metric;
  
  const double Avals[] = { 1, 2, 2, 5 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I( 1.0 );
  const MsqMatrix<2,2> A( Avals );
  const MsqMatrix<2,2> B( Bvals );
  
  MsqError err;
  MsqMatrix<2,2> g;
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
  
    // check inverted also for non-barier metrics
  if (Barrier) 
    return;
  
  const double Cvals[] = { -1.0, 0.5, 0.0, 1.0 };
  const MsqMatrix<2,2> C( Cvals );
  valid = metric.evaluate_with_grad( C, I, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate( C, I, v, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v, gv, 1e-6 );
}

template <class Metric> 
void Target2DTest<Metric>::compare_eval_with_grad_and_eval_with_hess()
{
  Metric metric;
  
  const double Avals[] = { 1, 2, 2, 5 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I( 1.0 );
  const MsqMatrix<2,2> A( Avals );
  const MsqMatrix<2,2> B( Bvals );
  
  MsqError err;
  MsqMatrix<2,2> g, h, hess[3];
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
  
    // check inverted also for non-barier metrics
  if (Barrier) 
    return;
  
  const double Cvals[] = { -1.0, 0.5, 0.0, 1.0 };
  const MsqMatrix<2,2> C( Cvals );
  valid = metric.evaluate_with_grad( C, I, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( C, I, hv, h, hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( gv, hv, 1e-6 );
  ASSERT_MATRICES_EQUAL( g, h, 1e-5 );
}

double eps( double mu ) { return std::max( mu*1e-2, 5e-2 ); }

template <class Metric> 
void Target2DTest<Metric>::compare_anaytic_and_numeric_grads()
{
  const double EPS_VAL = 1e-6;
  Metric metric;
  
  const double Avals[] = { 1, 2, 2, 5 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I( 1.0 );
  const MsqMatrix<2,2> A( Avals );
  const MsqMatrix<2,2> B( Bvals );
  
  MsqError err;
  MsqMatrix<2,2> num, ana;
  bool valid;
  double nval, aval;
  
  valid = metric.TargetMetric2D::evaluate_with_grad( I, A, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( I, A, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps(aval) );
  
  valid = metric.TargetMetric2D::evaluate_with_grad( A, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps(aval) );
  
  valid = metric.TargetMetric2D::evaluate_with_grad( I, B, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( I, B, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps(aval) );
  
  valid = metric.TargetMetric2D::evaluate_with_grad( B, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( B, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps(aval) );
   
  valid = metric.TargetMetric2D::evaluate_with_grad( A, B, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, B, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps(aval) );
  
  valid = metric.TargetMetric2D::evaluate_with_grad( B, A, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( B, A, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps(aval) );
  
    // check inverted also for non-barier metrics
  if (Barrier) 
    return;
  
  const double Cvals[] = { -1.0, 0.5, 0.0, 1.0 };
  const MsqMatrix<2,2> C( Cvals );
  valid = metric.TargetMetric2D::evaluate_with_grad( C, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( C, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps(aval) );
}

double epsh( double mu ) { return std::max( mu*1e-1, 1e-1 ); }

template <class Metric> 
void Target2DTest<Metric>::compare_anaytic_and_numeric_hess()
{
  const double EPS_VAL = 1e-6;
  Metric metric;
  
  const double Avals[] = { 1, 2, 2, 5 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I( 1.0 );
  const MsqMatrix<2,2> A( Avals );
  const MsqMatrix<2,2> B( Bvals );
  
  MsqError err;
  MsqMatrix<2,2> dmdA_num, dmdA_ana, d2mdA2_num[3], d2mdA2_ana[3];
  bool valid;
  double val_num, val_ana;
  
  valid = metric.TargetMetric2D::evaluate_with_hess( I, A, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, A, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], epsh(val_ana) );
  
  valid = metric.TargetMetric2D::evaluate_with_hess( A, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( A, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], epsh(val_ana) );
  
  valid = metric.TargetMetric2D::evaluate_with_hess( I, B, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( I, B, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], epsh(val_ana) );
  
  valid = metric.TargetMetric2D::evaluate_with_hess( B, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( B, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], epsh(val_ana) );
  
  valid = metric.TargetMetric2D::evaluate_with_hess( A, B, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( A, B, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], epsh(val_ana) );
  
  valid = metric.TargetMetric2D::evaluate_with_hess( B, A, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( B, A, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], epsh(val_ana) );
  
    // check inverted also for non-barier metrics
  if (Barrier) 
    return;
  
  const double Cvals[] = { -1.0, 0.5, 0.0, 1.0 };
  const MsqMatrix<2,2> C( Cvals );
  valid = metric.TargetMetric2D::evaluate_with_hess( C, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_hess( C, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], epsh(val_ana) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], epsh(val_ana) );
}


// implement metric that is the sum of the elements of 2A-W,
// such that the derivative of the result with resepct to
// each element of A is 2.
class GradTestMetric2D : public TargetMetric2D
{
  public:
    std::string get_name() const { return "GradTest"; }
  
    bool evaluate( const MsqMatrix<2,2>& A,
                   const MsqMatrix<2,2>& W,
                   double& result,
                   MsqError&  )
    {
      result = 0;
      for (int r = 0; r < 2; ++r) 
        for (int c = 0; c < 2; ++c)
          result += 2 * A(r,c) - W(r,c);
      return true;
    }
};

void TargetMetric2DTest::test_numerical_gradient()
{
  GradTestMetric2D metric;
  const double Avals[] = { 1, 2, 2, 5 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I( 1.0 );
  const MsqMatrix<2,2> A( Avals );
  const MsqMatrix<2,2> B( Bvals );
  
  MsqError err;
  MsqMatrix<2,2> d;
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
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  
  valid = metric.evaluate( A, I, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, I, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  
  valid = metric.evaluate( I, B, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( I, B, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  
  valid = metric.evaluate( B, I, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( B, I, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
   
  valid = metric.evaluate( A, B, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( A, B, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
  
  valid = metric.evaluate( B, A, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = metric.evaluate_with_grad( B, A, gval, d, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val, gval, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, d(1,1), 1e-6 );
}

// implement metric that is the |2A - A^T - W|^2
// such that the Hessian is the constant:
//  _          _
// | 2  0  0  0 |
// | 0 10 -8  0 |
// | 0 -8 10  0 |
// |_0  0  0  2_|

class HessTestMetric2D : public TargetMetric2D
{
  public:
    std::string get_name() const { return "HessTest"; }
  
    bool evaluate( const MsqMatrix<2,2>& A,
                   const MsqMatrix<2,2>& W,
                   double& result,
                   MsqError&  )
    {
      result = sqr_Frobenius(2*A - transpose(A) - W);
      return true;
    }
    bool evaluate_with_grad( const MsqMatrix<2,2>& A,
                             const MsqMatrix<2,2>& W,
                             double& result,
                             MsqMatrix<2,2>& wrt_A,
                             MsqError&  )
    {
      result = sqr_Frobenius(2*A - transpose(A) - W);
      wrt_A = 10*A - 8*transpose(A) - 4*W + 2*transpose(W);
      return true;
    }
};

void TargetMetric2DTest::test_numerical_hessian()
{
  HessTestMetric2D metric;
  const double Avals[] = { 1, 2, 2, 5 };
  const double Bvals[] = { -0.1, -0.15, -0.25, -0.8 };
  const MsqMatrix<2,2> I( 1.0 );
  const MsqMatrix<2,2> A( Avals );
  const MsqMatrix<2,2> B( Bvals );
  
  MsqError err;
  MsqMatrix<2,2> g,gh;
  MsqMatrix<2,2> h[3];
  bool valid;
  double val, hval;
  
  const double h_00[] = { 2, 0, 0, 10 };
  const double h_01[] = { 0, 0, -8, 0 };
  const double h_11[] = { 10, 0, 0, 2 };
  MsqMatrix<2,2> h00(h_00), h01(h_01), h11(h_11);
  
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
  ASSERT_MATRICES_EQUAL( h11, h[2], 1e-6 );
  
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
  ASSERT_MATRICES_EQUAL( h11, h[2], 1e-6 );
  
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
  ASSERT_MATRICES_EQUAL( h11, h[2], 1e-6 );
  
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
  ASSERT_MATRICES_EQUAL( h11, h[2], 1e-6 );
  
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
  ASSERT_MATRICES_EQUAL( h11, h[2], 1e-6 );
  
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
  ASSERT_MATRICES_EQUAL( h11, h[2], 1e-6 );
}



#include "TSquared2D.hpp"
#include "Target2DShape.hpp"
#include "Target2DShapeAlt1.hpp"
#include "Target2DShapeBarrier.hpp"
#include "Target2DShapeOrient.hpp"
#include "Target2DShapeOrientAlt1.hpp"
#include "Target2DShapeOrientBarrier.hpp"
#include "Target2DShapeOrientBarrierAlt1.hpp"
#include "Target2DShapeSize.hpp"
#include "Target2DShapeSizeAlt2.hpp"
#include "Target2DShapeSizeBarrier.hpp"
#include "Target2DShapeSizeBarrierAlt1.hpp"
#include "Target2DShapeSizeBarrierAlt2.hpp"
#include "Target2DShapeSizeOrient.hpp"
#include "Target2DShapeSizeOrientAlt1.hpp"
#include "Target2DShapeSizeOrientBarrier.hpp"
#include "Target2DShapeSizeOrientBarrierAlt1.hpp"
#include "InverseMeanRatio2D.hpp"
#include "Target2DSize.hpp"
#include "Target2DSizeBarrier.hpp"
#include "Target2DUntangle.hpp"
#include "Target2DUntangleAlt1.hpp"
#include "Target2DShapeSizeAlt1.hpp"

#define REGISTER_TARGET2D_TEST( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class METRIC ## Test : public Target2DTest<METRIC> { public: \
  METRIC ## Test () : Target2DTest<METRIC>( (SHAPE_INVAR), (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( METRIC ## Test ); \
  CPPUNIT_TEST (test_ideal_eval); \
  CPPUNIT_TEST (test_ideal_gradient); \
  CPPUNIT_TEST (test_inverted); \
  CPPUNIT_TEST (test_shape); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _FileRegister ("Target2DTest"); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _BaseRegister ( #METRIC "Test" )


// Macro arguments:
//  shape_invariant
//  size_invariant
//  orient_invariant
//  barrier
//  expected value for ideal element
#define REGISTER_TARGET2D_TEST_WITH_GRAD( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class METRIC ## Test : public Target2DTest<METRIC> { public: \
  METRIC ## Test () : Target2DTest<METRIC>( (SHAPE_INVAR), (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( METRIC ## Test ); \
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
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _FileRegister ("Target2DTest"); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _BaseRegister ( #METRIC "Test" )

#define REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class METRIC ## Test : public Target2DTest<METRIC> { public: \
  METRIC ## Test () : Target2DTest<METRIC>( (SHAPE_INVAR), (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( METRIC ## Test ); \
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
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _FileRegister ("Target2DTest"); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## Test > METRIC ## _BaseRegister ( #METRIC "Test" )

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric2DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric2DTest, "Target2DTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric2DTest, "TargetMetric2DTest" );

REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShape,                     false,  true,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeAlt1,                 false,  true,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeBarrier,              false,  true,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeOrient,               false,  true, false, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeOrientAlt1,           false,  true, false, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeOrientBarrier,        false,  true, false,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeOrientBarrierAlt1,    false,  true, false,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSize,                 false, false,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeAlt1,             false, false,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeAlt2,             false, false,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeBarrier,          false, false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeBarrierAlt1,      false, false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeBarrierAlt2,      false, false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeOrient,           false, false, false, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeOrientAlt1,       false, false, false, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeOrientBarrier,    false, false, false,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeOrientBarrierAlt1,false, false, false,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( InverseMeanRatio2D,                false,  true,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DSize,                       true, false,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DSizeBarrier,                true, false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DUntangle,                   true,  true,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DUntangleAlt1,               true,  true,  true, false, 0.0 );

class TSquared2DTest : public Target2DTest<TSquared2D> {
  public: 
    TSquared2DTest() : Target2DTest<TSquared2D>(false,false,false,false,0.0) {}
    CPPUNIT_TEST_SUITE( TSquared2DTest );
    CPPUNIT_TEST( compare_anaytic_and_numeric_grads );
    CPPUNIT_TEST( compare_eval_with_grad_and_eval_with_hess );
    CPPUNIT_TEST( compare_anaytic_and_numeric_hess );
    CPPUNIT_TEST_SUITE_END();
};
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, "Target2DTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, "TSquared2DTest" );
