/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

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

/*! \file JacobianMetricTest.cpp

Unit testing for the JacobianMetric class
\author Jasno Kraftcheck
*/


#include "JacobianMetric.hpp"
#include "TargetMetric2D.hpp"
#include "TargetMetric3D.hpp"
#include "UnitWeight.hpp"
#include "IdealTargetCalculator.hpp"
#include "MsqMatrix.hpp"
#include "QualityMetricTester.hpp"
#include "LinearFunctionSet.hpp"
#include "UnitUtil.hpp"
#include "PlanarDomain.hpp"
#include "SamplePoints.hpp"
#include "PatchData.hpp"

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else
#include <iostream>
#endif

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class JacobianMetricTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(JacobianMetricTest);
  CPPUNIT_TEST (test_negate_flag);
  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_evaluations);
  CPPUNIT_TEST (test_evaluate_2D);
  CPPUNIT_TEST (test_evaluate_3D);
  CPPUNIT_TEST (test_evaluate_2D_weight);
  CPPUNIT_TEST (test_evaluate_3D_weight);
  CPPUNIT_TEST (test_evaluate_2D_corner);
  CPPUNIT_TEST (test_evaluate_2D_edge);
  CPPUNIT_TEST (test_evaluate_2D_elem);
  CPPUNIT_TEST (test_evaluate_3D_corner);
  CPPUNIT_TEST (test_evaluate_3D_edge);
  CPPUNIT_TEST (test_evaluate_3D_edge);
  CPPUNIT_TEST (test_evaluate_3D_face);
  CPPUNIT_TEST (test_evaluate_with_indices);
  CPPUNIT_TEST (test_evaluate_fixed_indices);
  CPPUNIT_TEST_SUITE_END();
public:
  JacobianMetricTest() : 
    tester( QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON, &mf )
  {  }
  
  void test_negate_flag();
  void test_supported_types();
  void test_get_evaluations();
  void test_get_element_evaluations();
  void test_evaluate_2D();
  void test_evaluate_3D();
  void test_evaluate_2D_weight();
  void test_evaluate_3D_weight();
  void test_evaluate_2D_corner();
  void test_evaluate_2D_edge();
  void test_evaluate_2D_elem();
  void test_evaluate_3D_corner();
  void test_evaluate_3D_edge();
  void test_evaluate_3D_face();
  void test_evaluate_3D_elem();
  void test_evaluate_with_indices();
  void test_evaluate_fixed_indices();
  
private:
  void test_2d_eval_ortho_quad( unsigned dim );
  void test_3d_eval_ortho_hex( unsigned dim );

  LinearFunctionSet mf;
  QualityMetricTester tester;
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(JacobianMetricTest, "JacobianMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(JacobianMetricTest, "Unit");

template <typename B>
class FauxTarget : public B
{
public:
  int count;
  double value;
  bool rval;
  MsqMatrix<B::MATRIX_DIM,B::MATRIX_DIM> last_A, last_W;

  FauxTarget(double v) : count(0), value(v), rval(true) {}
  
  bool evaluate( const MsqMatrix<B::MATRIX_DIM,B::MATRIX_DIM>& A, 
                 const MsqMatrix<B::MATRIX_DIM,B::MATRIX_DIM>& W, 
                 double& result, MsqError&  )
  {
    last_A = A;
    last_W = W;
    result = value;
    ++count;
    return rval;
  }
};

class ScaleWeight : public WeightCalculator
{
  public:
    ScaleWeight( double s ) : value(s) {}
    double get_weight( PatchData&, size_t, const SamplePoints*, unsigned, MsqError& )
      { return value; }
    double value;
};

void JacobianMetricTest::test_negate_flag()
{
  SamplePoints pts( true, false, false, false );
  FauxTarget<TargetMetric2D> metric_2d(3.14159);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc;
  UnitWeight wc;
  JacobianMetric metric( &pts, &tc, &wc, &metric_2d, &metric_3d );
  
  CPPUNIT_ASSERT_EQUAL( 1, metric.get_negate_flag() );
}

void JacobianMetricTest::test_supported_types()
{
  SamplePoints pts(true);
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc;
  UnitWeight wc;
  
  JacobianMetric m( &pts, &tc, &wc, &metric_2d, &metric_3d );
  tester.test_supported_element_types( &m );
}

void JacobianMetricTest::test_get_evaluations()
{
  SamplePoints corners( true, false, false, false );
  SamplePoints corners_and_edges( true, true, false, false );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc;
  UnitWeight wc;
  
  JacobianMetric corner_metric( &corners, &tc, &wc, &metric_2d, &metric_3d );
  JacobianMetric edge_metric( &corners_and_edges, &tc, &wc, &metric_2d, &metric_3d );
  
  tester.test_get_sample_evaluations( &corner_metric, &corners );
  tester.test_get_sample_evaluations( &edge_metric, &corners_and_edges );
}
  

void JacobianMetricTest::test_get_element_evaluations()
{
  SamplePoints corners( true, false, false, false );
  SamplePoints edges( false, true, false, false );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc;
  UnitWeight wc;
  
  JacobianMetric corner_metric( &corners, &tc, &wc, &metric_2d, &metric_3d );
  JacobianMetric edge_metric( &edges, &tc, &wc, &metric_2d, &metric_3d );
  
  tester.test_get_in_element_evaluations( &corner_metric );
  tester.test_get_in_element_evaluations( &edge_metric );
}

static double col_dot_prod( MsqMatrix<2,2>& m )
  { return m(0,0) * m(0,1) + m(1,0) * m(1,1); }

void JacobianMetricTest::test_evaluate_2D()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  SamplePoints pts( true, false, false, false );
  FauxTarget<TargetMetric2D> metric_2d(3.14159);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc(false);
  UnitWeight wc;
  JacobianMetric m( &pts, &tc, &wc, &metric_2d, &metric_3d );
  
    // test with aligned elements
  metric_2d.count = metric_3d.count = 0;
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( metric_2d.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, metric_2d.count );
  CPPUNIT_ASSERT_EQUAL( 0, metric_3d.count );
  ASSERT_MATRICES_EQUAL( metric_2d.last_W, metric_2d.last_A, 1e-6 );
  
    // test that columns are orthogonal for ideal quad element
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, col_dot_prod(metric_2d.last_A), 1e-6 );
  
    // test with an element rotated about X-axis
  metric_2d.count = metric_3d.count = 0;
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  // rotate by 90 degrees about X axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    double z = pd.vertex_by_index(i)[2];
    pd.vertex_by_index(i)[2] = pd.vertex_by_index(i)[1];
    pd.vertex_by_index(i)[1] =-z;
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( metric_2d.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, metric_2d.count );
  CPPUNIT_ASSERT_EQUAL( 0, metric_3d.count );
  ASSERT_MATRICES_EQUAL( metric_2d.last_W, metric_2d.last_A, 1e-6 );
  
    // test with an element rotated about Y-axis
  metric_2d.count = metric_3d.count = 0;
  tester.get_ideal_element( TRIANGLE, true, pd );
  // rotate by -90 degrees about Y axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    double z = pd.vertex_by_index(i)[2];
    pd.vertex_by_index(i)[2] = -pd.vertex_by_index(i)[0];
    pd.vertex_by_index(i)[0] = z;
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( metric_2d.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, metric_2d.count );
  CPPUNIT_ASSERT_EQUAL( 0, metric_3d.count );
  ASSERT_MATRICES_EQUAL( metric_2d.last_W, metric_2d.last_A, 1e-6 );
}  
 
  
void JacobianMetricTest::test_evaluate_3D()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  SamplePoints pts( true, false, false, false );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(exp(1.0));
  IdealTargetCalculator tc(false);
  UnitWeight wc;
  JacobianMetric m( &pts, &tc, &wc, &metric_2d, &metric_3d );
  
    // test with aligned elements
  metric_2d.count = metric_3d.count = 0;
  tester.get_ideal_element( HEXAHEDRON, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( metric_3d.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 0, metric_2d.count );
  CPPUNIT_ASSERT_EQUAL( 1, metric_3d.count );
  ASSERT_MATRICES_EQUAL( metric_3d.last_W, metric_3d.last_A, 1e-6 );
  
    // test that columns are orthogonal for ideal hex element
  MsqMatrix<3,3> A = metric_3d.last_A;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(1) % A.column(2), 1e-6 );
  
    // test with rotated element
  metric_2d.count = metric_3d.count = 0;
  tester.get_ideal_element( TETRAHEDRON, true, pd );
  // rotate by 90-degrees about X axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    double z = pd.vertex_by_index(i)[2];
    pd.vertex_by_index(i)[2] = pd.vertex_by_index(i)[1];
    pd.vertex_by_index(i)[1] =-z;
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( metric_3d.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 0, metric_2d.count );
  CPPUNIT_ASSERT_EQUAL( 1, metric_3d.count );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( determinant(metric_3d.last_W),
                                determinant(metric_3d.last_A),
                                1e-6 );
}


void JacobianMetricTest::test_evaluate_2D_weight()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  SamplePoints pts( true, false, false, false );
  FauxTarget<TargetMetric2D> metric_2d(2.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc(false);
  ScaleWeight wc(3.14159);
  JacobianMetric m( &pts, &tc, &wc, &metric_2d, &metric_3d );
  
  tester.get_ideal_element( TRIANGLE, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( metric_2d.value*wc.value, value, DBL_EPSILON );
}


void JacobianMetricTest::test_evaluate_3D_weight()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  SamplePoints pts( true, false, false, false );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.5);
  IdealTargetCalculator tc(false);
  ScaleWeight wc(3.14159);
  JacobianMetric m( &pts, &tc, &wc, &metric_2d, &metric_3d );
  
  tester.get_ideal_element( PRISM, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( metric_3d.value*wc.value, value, DBL_EPSILON );
}

void JacobianMetricTest::test_2d_eval_ortho_quad( unsigned dim )
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  bool locs[4] = {false, false, false, false};
  locs[dim] = true;
  
  SamplePoints pts( locs[0], locs[1], locs[2], locs[3] );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc(false);
  UnitWeight wc;
  JacobianMetric m( &pts, &tc, &wc, &metric_2d, &metric_3d );
  
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( 1, metric_2d.count );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, col_dot_prod(metric_2d.last_A), DBL_EPSILON );
}  

void JacobianMetricTest::test_3d_eval_ortho_hex( unsigned dim )
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  bool locs[4] = {false, false, false, false};
  locs[dim] = true;
  
  SamplePoints pts( locs[0], locs[1], locs[2], locs[3] );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc(false);
  UnitWeight wc;
  JacobianMetric m( &pts, &tc, &wc, &metric_2d, &metric_3d );
  
  tester.get_ideal_element( HEXAHEDRON, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( 1, metric_3d.count );
  
    // test that columns are orthogonal for ideal hex element
  MsqMatrix<3,3> A = metric_3d.last_A;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(1) % A.column(2), 1e-6 );
}  

void JacobianMetricTest::test_evaluate_2D_corner()
  { test_2d_eval_ortho_quad(0); }

void JacobianMetricTest::test_evaluate_2D_edge()
  { test_2d_eval_ortho_quad(1); }

void JacobianMetricTest::test_evaluate_2D_elem()
  { test_2d_eval_ortho_quad(3); }

void JacobianMetricTest::test_evaluate_3D_corner()
  { test_3d_eval_ortho_hex(0); }
  
void JacobianMetricTest::test_evaluate_3D_edge()
  { test_3d_eval_ortho_hex(1); }
  
void JacobianMetricTest::test_evaluate_3D_face()
  { test_3d_eval_ortho_hex(2); }
  
void JacobianMetricTest::test_evaluate_3D_elem()
  { test_3d_eval_ortho_hex(3); }

void JacobianMetricTest::test_evaluate_with_indices()
{
  SamplePoints corners( true, false, false, false );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc;
  UnitWeight wc;
  
  JacobianMetric m( &corners, &tc, &wc, &metric_2d, &metric_3d );
  tester.test_get_sample_indices( &m, &corners );
  tester.compare_eval_and_eval_with_indices( &m );
}

void JacobianMetricTest::test_evaluate_fixed_indices()
{
  SamplePoints corners( true, false, false, false );
  FauxTarget<TargetMetric2D> metric_2d(0.0);
  FauxTarget<TargetMetric3D> metric_3d(0.0);
  IdealTargetCalculator tc;
  UnitWeight wc;
  
  JacobianMetric m( &corners, &tc, &wc, &metric_2d, &metric_3d );
  tester.test_get_indices_fixed( &m );
}

