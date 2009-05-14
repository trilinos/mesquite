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
#include "QualityMetricTester.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "TargetMetric2D.hpp"
#include "TMPQualityMetric.hpp"
#include "IdealTargetCalculator.hpp"
#include "Settings.hpp"
#include "UnitUtil.hpp"

static const EntityTopology SurfElems[] = { TRIANGLE, QUADRILATERAL };

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
  Settings settings;
  QualityMetricTester tester;
  IdealTargetCalculator target;
  Metric test_metric;
  TMPQualityMetric metric;
  bool sizeInvariant, orientInvariant, Barrier;
  double idealVal;
public:
  Target2DTest( bool size_invariant, bool orient_invariant, bool barrier, double ideal_element_val )
    : tester( SurfElems, sizeof(SurfElems)/sizeof(SurfElems[0]), &settings ),
      metric( &target, &test_metric, 0 ),
      sizeInvariant(size_invariant), orientInvariant(orient_invariant), Barrier(barrier),
      idealVal(ideal_element_val)
    {}
  
  inline void test_ideal_element_eval() {
    tester.test_evaluate_unit_element( &metric, TRIANGLE, idealVal );
    tester.test_evaluate_unit_element( &metric, QUADRILATERAL, idealVal );
  }
  
  inline void test_ideal_element_gradient() {
    tester.test_ideal_element_zero_gradient( &metric, true );
  }

  inline void test_inverted_element_eval() {
    tester.test_evaluate_inverted_element( &metric, !Barrier );
  }
  
  inline void test_measures_quality() {
    if (sizeInvariant && orientInvariant)
      tester.test_measures_quality( &metric );
  }
  
  inline void test_location_invariant() {
    tester.test_location_invariant( &metric);
    tester.test_grad_location_invariant( &metric );
  }
  
  inline void test_scale() {
    if (sizeInvariant) {
      // these tests is not applicable to the target metrics.
      //tester.test_scale_invariant( &metric );
    }
    else {
      tester.test_measures_size( &metric, true );
    }
  }
  
  inline void test_orient() {
    if (orientInvariant) {
      tester.test_orient_invariant( &metric );
      tester.test_grad_orient_invariant( &metric );
    }
    else {
      tester.test_measures_in_plane_orientation( &metric );
    }
  }
  
  void compare_anaytic_and_numeric_grads();
  void compare_anaytic_and_numeric_hess();
  void compare_eval_and_eval_with_grad();
  void compare_eval_with_grad_and_eval_with_hess();
};

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
}

double eps( double mu ) { return std::max( mu*1e-3, 5e-3 ); }

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
}


// implement metric that is the sum of the elements of 2A-W,
// such that the derivative of the result with resepct to
// each element of A is 2.
class GradTestMetric2D : public TargetMetric2D
{
  public:
    msq_std::string get_name() const { return "GradTest"; }
  
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
    msq_std::string get_name() const { return "HessTest"; }
  
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
#include "Target2DShapeBarrier.hpp"
#include "Target2DShapeOrient.hpp"
#include "Target2DShapeOrientAlt1.hpp"
#include "Target2DShapeOrientAlt2.hpp"
#include "Target2DShapeOrientBarrier.hpp"
#include "Target2DShapeSize.hpp"
#include "Target2DShapeSizeBarrier.hpp"
#include "Target2DShapeSizeBarrierAlt1.hpp"
#include "Target2DShapeSizeBarrierAlt2.hpp"
#include "Target2DShapeSizeOrient.hpp"
#include "Target2DShapeSizeOrientAlt1.hpp"
#include "Target2DShapeSizeOrientBarrier.hpp"
#include "Target2DShapeSizeOrientBarrierAlt2.hpp"
#include "InverseMeanRatio2D.hpp"
#include "Target2DSize.hpp"
#include "Target2DSizeBarrier.hpp"

#define REGISTER_TARGET2D_TEST( METRIC, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class Test_ ## METRIC : public Target2DTest<METRIC> { public: \
  Test_ ## METRIC () : Target2DTest<METRIC>( (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( Test_ ## METRIC ); \
  CPPUNIT_TEST (test_ideal_element_eval); \
  CPPUNIT_TEST (test_ideal_element_gradient); \
  CPPUNIT_TEST (test_inverted_element_eval); \
  CPPUNIT_TEST (test_measures_quality); \
  CPPUNIT_TEST (test_location_invariant); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _FileRegister ("Target2DTest"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _BaseRegister ( "Test_" #METRIC )


// Macro arguments:
//  size_invarient
//  orient_invarient
//  barrier
//  expected value for ideal element
#define REGISTER_TARGET2D_TEST_WITH_GRAD( METRIC, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class Test_ ## METRIC : public Target2DTest<METRIC> { public: \
  Test_ ## METRIC () : Target2DTest<METRIC>( (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( Test_ ## METRIC ); \
  CPPUNIT_TEST (test_ideal_element_eval); \
  CPPUNIT_TEST (test_ideal_element_gradient); \
  CPPUNIT_TEST (test_inverted_element_eval); \
  CPPUNIT_TEST (test_measures_quality); \
  CPPUNIT_TEST (test_location_invariant); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_grads); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _FileRegister ("Target2DTest"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _BaseRegister ( "Test_" #METRIC )

#define REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( METRIC, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL_VAL ) \
class Test_ ## METRIC : public Target2DTest<METRIC> { public: \
  Test_ ## METRIC () : Target2DTest<METRIC>( (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), (IDEAL_VAL) ) {} \
  CPPUNIT_TEST_SUITE( Test_ ## METRIC ); \
  CPPUNIT_TEST (test_ideal_element_eval); \
  CPPUNIT_TEST (test_ideal_element_gradient); \
  CPPUNIT_TEST (test_inverted_element_eval); \
  CPPUNIT_TEST (test_measures_quality); \
  CPPUNIT_TEST (test_location_invariant); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_grads); \
  CPPUNIT_TEST (compare_eval_with_grad_and_eval_with_hess); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_hess); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _FileRegister ("Target2DTest"); \
CPPUNIT_NS::AutoRegisterSuite< Test_ ## METRIC > METRIC ## _BaseRegister ( "Test_" #METRIC )

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric2DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric2DTest, "Target2DTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TargetMetric2DTest, "TargetMetric2DTest" );

REGISTER_TARGET2D_TEST               ( Target2DShape,                      true,  true, false, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeBarrier,               true,  true,  true, 1.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeOrient,                true, false, false, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeOrientAlt1,            true, false, false, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeOrientAlt2,            true, false, false, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeOrientBarrier,         true, false,  true, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeSize,                 false,  true, false, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeSizeBarrier,          false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeSizeBarrierAlt1,      false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeSizeBarrierAlt2,      false,  true,  true, 1.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DShapeSizeOrient,           false, false, false, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeSizeOrientAlt1,       false, false, false, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeSizeOrientBarrier,    false, false,  true, 0.0 );
REGISTER_TARGET2D_TEST               ( Target2DShapeSizeOrientBarrierAlt2,false, false,  true, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( InverseMeanRatio2D,                 true,  true,  true, 1.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DSize,                      false,  true, false, 0.0 );
REGISTER_TARGET2D_TEST_WITH_2ND_DERIV( Target2DSizeBarrier,               false,  true,  true, 2.0 );

class Test_TSquared2D : public Target2DTest<TSquared2D> {
  public: 
    Test_TSquared2D() : Target2DTest<TSquared2D>(false,false,false,0.0) {}
    CPPUNIT_TEST_SUITE( Test_TSquared2D );
    CPPUNIT_TEST( compare_anaytic_and_numeric_grads );
    CPPUNIT_TEST( compare_eval_with_grad_and_eval_with_hess );
    CPPUNIT_TEST( compare_anaytic_and_numeric_hess );
    CPPUNIT_TEST_SUITE_END();
};
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( Test_TSquared2D, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( Test_TSquared2D, "Target2DTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( Test_TSquared2D, "Test_TSquared2D" );
