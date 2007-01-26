/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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

  ***************************************************************** */


/** \file PowerMeanPTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "PowerMeanP.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "UnitUtil.hpp"
#include "VertexQM.hpp"
#include "MsqHessian.hpp"

#include "ObjectiveFunctionTests.hpp"


using namespace Mesquite;
using namespace std;

const double EPSILON = 1e-4;

class PowerMeanPTest : public CppUnit::TestFixture
{
private:
  void test_eval_type( OFTestMode eval_func, ObjectiveFunction::EvalType type );
  void test_evaluate( double power );
  void test_gradient( double power );
  void test_Hessian( double power );
  
  void check_result( PatchData& pd, double power, double value, 
                     Vector3D* gradient = 0, Matrix3D* Hessian = 0 );

  CPPUNIT_TEST_SUITE( PowerMeanPTest );

  CPPUNIT_TEST( test_eval_calc );
  CPPUNIT_TEST( test_eval_accum );
  CPPUNIT_TEST( test_eval_save );
  CPPUNIT_TEST( test_eval_update );
  CPPUNIT_TEST( test_eval_temp );

  CPPUNIT_TEST( test_grad_calc );
  CPPUNIT_TEST( test_grad_save );
  CPPUNIT_TEST( test_grad_update );
  CPPUNIT_TEST( test_grad_temp );

  CPPUNIT_TEST( test_Hess_calc );
  CPPUNIT_TEST( test_Hess_save );
  CPPUNIT_TEST( test_Hess_update );
  CPPUNIT_TEST( test_Hess_temp );
  
  CPPUNIT_TEST( test_clone );
  
  CPPUNIT_TEST( test_failed_metric_in_eval );
  CPPUNIT_TEST( test_failed_metric_in_grad );
  CPPUNIT_TEST( test_failed_metric_in_Hess );
  
  CPPUNIT_TEST( test_false_metric_in_eval );
  CPPUNIT_TEST( test_false_metric_in_grad );
  CPPUNIT_TEST( test_false_metric_in_Hess );
  
  CPPUNIT_TEST( test_evaluate_arithmatic );
  CPPUNIT_TEST( test_evaluate_rms );
  
  CPPUNIT_TEST( test_gradient_arithmatic );
  CPPUNIT_TEST( test_gradient_rms );

  CPPUNIT_TEST( test_Hessian_arithmatic );
  CPPUNIT_TEST( test_Hessian_rms );
  
  CPPUNIT_TEST( compare_gradient_arithmatic );
  CPPUNIT_TEST( compare_gradient_rms );
 
  CPPUNIT_TEST( compare_hessian_gradient_arithmatic );
  CPPUNIT_TEST( compare_hessian_gradient_rms );
  
  CPPUNIT_TEST( test_negate_eval );
  CPPUNIT_TEST( test_negate_grad );
  CPPUNIT_TEST( test_negate_hess );

  CPPUNIT_TEST_SUITE_END();
  PatchData mPatch;

public:

  void setUp();

  void test_eval_calc()   { test_eval_type( EVAL, ObjectiveFunction::CALCULATE ); }
  void test_eval_accum()  { test_eval_type( EVAL, ObjectiveFunction::ACCUMULATE ); }
  void test_eval_save()   { test_eval_type( EVAL, ObjectiveFunction::SAVE ); }
  void test_eval_update() { test_eval_type( EVAL, ObjectiveFunction::UPDATE ); }
  void test_eval_temp()   { test_eval_type( EVAL, ObjectiveFunction::TEMPORARY ); }

  void test_grad_calc()   { test_eval_type( GRAD, ObjectiveFunction::CALCULATE ); }
  void test_grad_save()   { test_eval_type( GRAD, ObjectiveFunction::SAVE ); }
  void test_grad_update() { test_eval_type( GRAD, ObjectiveFunction::UPDATE ); }
  void test_grad_temp()   { test_eval_type( GRAD, ObjectiveFunction::TEMPORARY ); }

  void test_Hess_calc()   { test_eval_type( HESS, ObjectiveFunction::CALCULATE ); }
  void test_Hess_save()   { test_eval_type( HESS, ObjectiveFunction::SAVE ); }
  void test_Hess_update() { test_eval_type( HESS, ObjectiveFunction::UPDATE ); }
  void test_Hess_temp()   { test_eval_type( HESS, ObjectiveFunction::TEMPORARY ); }

  void test_clone() { PowerMeanP of( 1, NULL ); ::test_clone(&of); }
  
  void test_failed_metric_in_eval() 
    { PowerMeanP of( 1, NULL ); test_handles_qm_error( EVAL, &of); }
  void test_failed_metric_in_grad() 
    { PowerMeanP of( 1, NULL ); test_handles_qm_error( GRAD, &of); }
  void test_failed_metric_in_Hess() 
    { PowerMeanP of( 1, NULL ); test_handles_qm_error( HESS, &of); }
  
  void test_false_metric_in_eval() 
    { PowerMeanP of( 1, NULL ); test_handles_invalid_qm( EVAL, &of); }
  void test_false_metric_in_grad() 
    { PowerMeanP of( 1, NULL ); test_handles_invalid_qm( GRAD, &of); }
  void test_false_metric_in_Hess() 
    { PowerMeanP of( 1, NULL ); test_handles_invalid_qm( HESS, &of); }
  
  void test_evaluate_arithmatic() { test_evaluate( 1 ); }
  void test_evaluate_rms()        { test_evaluate( 2 ); }
  
  void test_gradient_arithmatic() { test_gradient( 1 ); }
  void test_gradient_rms()        { test_gradient( 2 ); }

  void test_Hessian_arithmatic()  { test_Hessian( 1 ); }
  void test_Hessian_rms()         { test_Hessian( 2 ); }
  
  void compare_gradient_arithmatic() 
    { PowerMeanP of( 1, NULL ); compare_numerical_gradient( &of ); }
  void compare_gradient_rms()
    { PowerMeanP of( 2, NULL ); compare_numerical_gradient( &of ); }
  
  void compare_hessian_gradient_arithmatic() 
    { PowerMeanP of( 1, NULL ); compare_hessian_gradient( &of ); }
  void compare_hessian_gradient_rms()
    { PowerMeanP of( 2, NULL ); compare_hessian_gradient( &of ); }
    
  void test_negate_eval()
    { PowerMeanP of( 2, NULL ); test_negate_flag( EVAL, &of ); }
  void test_negate_grad()
    { PowerMeanP of( 2, NULL ); test_negate_flag( GRAD, &of ); }
  void test_negate_hess()
    { PowerMeanP of( 2, NULL ); test_negate_flag( HESS, &of ); }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerMeanPTest, "PowerMeanPTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerMeanPTest, "Unit");

void PowerMeanPTest::setUp()
{
  MsqPrintError err(std::cout);
  
  // Create a triangle mesh with three free vertices
  const double coords[] = { 1, 1, 0,
                            2, 1, 0,
                            3, 1, 0,
                            2, 2, 0,
                            1, 3, 0,
                            1, 2, 0 };
  const bool fixed_vtx[] = { true,
                             false,
                             true,
                             false,
                             true,
                             false };
  const size_t tri_conn[] = { 0, 1, 5,
                              1, 2, 3,
                              3, 4, 5,
                              1, 3, 5 };
  mPatch.fill( 6, coords, 
               4, TRIANGLE, tri_conn,
               fixed_vtx, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));

}

/** Define a fake quality metric for testing the objective function
 *
 * Returns, for each vertex in the patch, the distance of that
 * vertex from the origin as the quality.  Each evaluation depends
 * only on a single vertex.
 */
class DistTestMetric : public VertexQM
{
public:
  DistTestMetric() : falseEval(false), failEval(false) {}
  string get_name() const { return "Fake metric for testing objective function"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t vtx_idx, double &value, MsqError& err );
  bool evaluate_with_indices( PatchData& pd, size_t vtx_idx, double &value, vector<size_t>& indices, MsqError& err );
  //bool evaluate_with_gradient( PatchData& pd, size_t vtx_idx, double &value, vector<size_t>& indices, vector<Vector3D>& grad, MsqError& err );
  bool falseEval;
  bool failEval;
};

bool DistTestMetric::evaluate( PatchData& pd, size_t vtx_idx, 
                               double &value, MsqError& err )
{
  if (failEval) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return true;
  }

  MsqVertex& vtx = pd.vertex_by_index( vtx_idx );
  value = vtx.length_squared();
  return !falseEval;
}

bool DistTestMetric::evaluate_with_indices( PatchData& pd, size_t vtx_idx, 
                               double &value, vector<size_t>& indices, 
                               MsqError& err )
{
  if (failEval) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return true;
  }

  indices.clear();
  if (vtx_idx < pd.num_free_vertices())
    indices.push_back( vtx_idx );
  
  MsqVertex& vtx = pd.vertex_by_index( vtx_idx );
  value = vtx.length_squared();
  return !falseEval;
}


void PowerMeanPTest::test_eval_type( OFTestMode eval_func, ObjectiveFunction::EvalType type )
{
  PowerMeanP func( 1, NULL );
  ::test_eval_type( type, eval_func, &func );
}


void PowerMeanPTest::test_evaluate( double power )
{
  MsqPrintError err(cout);
  double value;
  bool rval;
  
  DistTestMetric metric;
  PowerMeanP func( power, &metric );
  rval = func.evaluate( ObjectiveFunction::CALCULATE, mPatch, value, OF_FREE_EVALS_ONLY, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  
  check_result( mPatch, power, value );
}

void PowerMeanPTest::test_gradient( double power )
{
  MsqPrintError err(cout);
  double value;
  bool rval;
  vector<Vector3D> grad;
  
  DistTestMetric metric;
  PowerMeanP func( power, &metric );
  rval = func.evaluate_with_gradient( ObjectiveFunction::CALCULATE, mPatch, value, grad, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL(mPatch.num_free_vertices(), grad.size());
  
  check_result( mPatch, power, value, &grad[0] );
}


void PowerMeanPTest::test_Hessian( double power )
{
  MsqPrintError err(cout);
  double value;
  bool rval;
  vector<Vector3D> grad;
  MsqHessian Hess;
  Hess.initialize( mPatch, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  
  DistTestMetric metric;
  PowerMeanP func( power, &metric );
  rval = func.evaluate_with_Hessian( ObjectiveFunction::CALCULATE, mPatch, value, grad, Hess, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  size_t n = mPatch.num_free_vertices();
  CPPUNIT_ASSERT_EQUAL( n, grad.size() );
  CPPUNIT_ASSERT_EQUAL( n, Hess.size() );
  
  Matrix3D zero( 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  vector<Matrix3D> Hessians(n);
  for (size_t r = 0; r < n; ++r) {
    Matrix3D* mat = Hess.get_block( r, r );
    CPPUNIT_ASSERT( mat != 0 );
    Hessians[r] = *mat;
    
    for (size_t c = r+1; c < n; ++c) {
      mat = Hess.get_block( r, c );
      if (mat)
        CPPUNIT_ASSERT_MATRICES_EQUAL( zero, *mat, EPSILON );
    }
  }
 
  check_result( mPatch, power, value, &grad[0], &Hessians[0] );
}

void PowerMeanPTest::check_result( PatchData& pd, double power, double value, 
                                   Vector3D* gradient, Matrix3D* Hessian )
{
  MsqPrintError err(cout);
  double mvalue, sum = 0;
  bool rval;
  vector<Vector3D> grads;
  vector<Matrix3D> Hess;
  vector<size_t> indices;
  
  
  DistTestMetric metric;
  
  size_t N = pd.num_free_vertices();
  for (size_t i = 0; i < N; ++i)
  {
    rval = metric.evaluate_with_Hessian( pd, i, mvalue, indices, grads, Hess, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) && rval );
    sum += pow( mvalue, power );
    
    if (!OF_FREE_EVALS_ONLY && indices.empty())
      continue;
      
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() );
    CPPUNIT_ASSERT_EQUAL( i, indices[0] );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, grads.size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, Hess.size() );
    
    if (gradient)
    {
      double f = power * pow( mvalue, power - 1 ) / N;
      CPPUNIT_ASSERT_VECTORS_EQUAL( f * grads[0], gradient[i], EPSILON );
    }
    if (Hessian)
    {
      double f = power / N;
      double p2 = (power - 1) * pow( mvalue, power - 2 );
      double p1 = pow( mvalue, power - 1 );
      Matrix3D m;
      m.outer_product( grads[0], grads[0] );
      m *= p2;
      m += p1 * Hess[0];
      m *= f;
      CPPUNIT_ASSERT_MATRICES_EQUAL( m, Hessian[i], EPSILON );
    }  
  }
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( sum / N, value, EPSILON );
}

