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

/*! \file TMPQualityMetricTest.cpp

Unit testing for the TMPQualityMetric class
\author Jasno Kraftcheck
*/


#include "TMPQualityMetric.hpp"
#include "TargetMetric2D.hpp"
#include "TargetMetric3D.hpp"
#include "IdealTargetCalculator.hpp"
#include "MsqMatrix.hpp"
#include "QualityMetricTester.hpp"
#include "LinearFunctionSet.hpp"
#include "UnitUtil.hpp"
#include "PlanarDomain.hpp"
#include "SamplePoints.hpp"
#include "PatchData.hpp"
#include "Target3DShapeSizeOrient.hpp"
#include "Target2DShapeSizeOrient.hpp"
#include "WeightCalculator.hpp"
#include "InverseMeanRatio2D.hpp"
#include "InverseMeanRatio3D.hpp"
#include "ElementPMeanP.hpp"
#include "IdealWeightInverseMeanRatio.hpp"

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else
#include <iostream>
#endif

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;


/** Target metric (templatized by dimension) for use in misc. tests.
 *  'evaluate' method records input values and returns a constant.
 */
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

/** Weight calculator used for testing.  Returns constant weight. */
class ScaleWeight : public WeightCalculator
{
  public:
    ScaleWeight( double s ) : value(s) {}
    double get_weight( PatchData&, size_t, const SamplePoints*, unsigned, MsqError& )
      { return value; }
    double value;
};

/** wrapper class to force numeric approximation of derivatives */
template <class Base>
class NumericalTarget : public Base
{
public:
  
  NumericalTarget( Base* real_metric ) : mMetric(real_metric) {}

  ~NumericalTarget() {}

  bool evaluate( const MsqMatrix<Base::MATRIX_DIM,Base::MATRIX_DIM>& A, 
                 const MsqMatrix<Base::MATRIX_DIM,Base::MATRIX_DIM>& W, 
                 double& result, 
                 MsqError& err )
  { return mMetric->evaluate( A, W, result, err ); }
private:
  Base* mMetric;
};


/** Simple target metric for testing first partial derivatives.  
 *  \f$\mu(A,W) = |A|^2\f$
 *  \f$\frac{\partial\mu}{\partial \A} = 2 A \f$
 */
class TestGradTargetMetric3D : public TargetMetric3D
{
  public:
  
    bool evaluate( const MsqMatrix<3,3>& A, const MsqMatrix<3,3>&, double& result, MsqError& err )
      { result = sqr_Frobenius(A); return true; }
    
    bool evaluate_with_grad( const MsqMatrix<3,3>& A, 
                             const MsqMatrix<3,3>&,
                             double& result,
                             MsqMatrix<3,3>& d,
                             MsqError& err )
    {
      result = sqr_Frobenius(A);
      d = 2*A;
      return true;
    }
};
class TestGradTargetMetric2D : public TargetMetric2D
{
  public:
  
    bool evaluate( const MsqMatrix<2,2>& A, const MsqMatrix<2,2>&, double& result, MsqError& err )
      { result = sqr_Frobenius(A); return true; }
    
    bool evaluate_with_grad( const MsqMatrix<2,2>& A, 
                             const MsqMatrix<2,2>&,
                             double& result,
                             MsqMatrix<2,2>& d,
                             MsqError& err )
    {
      result = sqr_Frobenius(A);
      d = 2*A;
      return true;
    }
};

class TMPQualityMetricTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(TMPQualityMetricTest);
  
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
  
  CPPUNIT_TEST (test_sample_indices);
  CPPUNIT_TEST (test_evaluate_with_indices);
  CPPUNIT_TEST (test_evaluate_fixed_indices);
  
  CPPUNIT_TEST (test_gradient_2D);
  CPPUNIT_TEST (test_gradient_3D);
  CPPUNIT_TEST (compare_indices_and_gradient);
  CPPUNIT_TEST (test_ideal_element_gradient);
  CPPUNIT_TEST (compare_analytical_and_numerical_gradient);
  CPPUNIT_TEST (test_weighted_gradients);
  CPPUNIT_TEST (test_gradient_with_fixed_vertices);

  CPPUNIT_TEST (compare_indices_and_hessian);
  CPPUNIT_TEST (compare_gradient_and_hessian);
  CPPUNIT_TEST (compare_analytical_and_numerical_hessians);
  CPPUNIT_TEST (test_symmetric_hessian_diagonal);
  CPPUNIT_TEST (test_weighted_hessians);
  CPPUNIT_TEST (test_hessian_with_fixed_vertices);

  CPPUNIT_TEST (compare_indices_and_diagonal);
  CPPUNIT_TEST (compare_gradient_and_diagonal);
  CPPUNIT_TEST (compare_analytical_and_numerical_diagonals);
  CPPUNIT_TEST (test_weighted_diagonals);
  CPPUNIT_TEST (test_diagonal_with_fixed_vertices);
  
  CPPUNIT_TEST (test_inverse_mean_ratio_grad);
  CPPUNIT_TEST (test_inverse_mean_ratio_hess);
  CPPUNIT_TEST (test_inverse_mean_ratio_hess_diag);
  CPPUNIT_TEST (regression_inverse_mean_ratio_grad);
  CPPUNIT_TEST (regression_inverse_mean_ratio_hess);

  CPPUNIT_TEST_SUITE_END();
  
  void test_2d_eval_ortho_quad( unsigned dim );
  void test_3d_eval_ortho_hex( unsigned dim );

  QualityMetricTester tester;

  LinearFunctionSet mf;
  SamplePoints corners, center;
  IdealTargetCalculator ideal;
  ScaleWeight e_weight;

  FauxTarget<TargetMetric2D> faux_2d_pi, faux_2d_zero;
  FauxTarget<TargetMetric3D> faux_3d_zero, faux_3d_two;
  Target3DShapeSizeOrient test_metric_3D;
  Target2DShapeSizeOrient test_metric_2D;
  NumericalTarget<TargetMetric3D> num_metric_3D;
  NumericalTarget<TargetMetric2D> num_metric_2D;
  TMPQualityMetric test_qm, zero_qm, weight_qm, center_qm;
  
public:
  TMPQualityMetricTest() : 
    tester( QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON, &mf ),
    corners( true, false, false, false ),
    center( false, false, false, true ),
    e_weight( 2.7182818284590451 ),
    faux_2d_pi(3.14159), faux_2d_zero(0.0),
    faux_3d_zero(0.0), faux_3d_two(2.0),
    num_metric_3D( &test_metric_3D ),
    num_metric_2D( &test_metric_2D ),
    test_qm( &corners, &ideal, &num_metric_2D, &num_metric_3D ),
    zero_qm( &corners, &ideal, &faux_2d_zero, &faux_3d_zero ),
    weight_qm( &corners, &ideal, &e_weight, &test_metric_2D, &test_metric_3D ),
    center_qm( &center, &ideal, &test_metric_2D, &test_metric_3D )
  {  }
  
  void test_negate_flag()
    { CPPUNIT_ASSERT_EQUAL( 1, zero_qm.get_negate_flag() ); }
  void test_supported_types()     
    { tester.test_supported_element_types( &zero_qm ); }
  void test_get_evaluations();
  void test_get_element_evaluations();
  void test_evaluate_2D();
  void test_evaluate_3D();
  void test_evaluate_2D_weight();
  void test_evaluate_3D_weight();
  void test_evaluate_2D_corner()  { test_2d_eval_ortho_quad(0); }
  void test_evaluate_2D_edge()    { test_2d_eval_ortho_quad(1); }
  void test_evaluate_2D_elem()    { test_2d_eval_ortho_quad(2); }
  void test_evaluate_3D_corner()  { test_3d_eval_ortho_hex(0); }
  void test_evaluate_3D_edge()    { test_3d_eval_ortho_hex(1); }
  void test_evaluate_3D_face()    { test_3d_eval_ortho_hex(2); }
  void test_evaluate_3D_elem()    { test_3d_eval_ortho_hex(3); }
  void test_gradient_2D();
  void test_gradient_3D();

  void test_sample_indices()
    { tester.test_get_sample_indices( &zero_qm, &corners ); }
  void test_evaluate_with_indices()
    { tester.compare_eval_and_eval_with_indices( &zero_qm ); }
  void test_evaluate_fixed_indices() 
    { tester.test_get_indices_fixed( &zero_qm ); }
    
  void compare_indices_and_gradient()
    { tester.compare_eval_with_indices_and_eval_with_gradient( &test_qm ); }
  void test_ideal_element_gradient()
    { tester.test_ideal_element_zero_gradient( &test_qm, true ); }
  void compare_analytical_and_numerical_gradient()
    { compare_analytical_and_numerical_gradients( &test_qm ); }
  void test_weighted_gradients()
    { compare_analytical_and_numerical_gradients( &weight_qm ); }
  void test_gradient_with_fixed_vertices()
    { tester.test_gradient_with_fixed_vertex( &center_qm ); }

  void compare_indices_and_hessian()
    { tester.compare_eval_with_indices_and_eval_with_hessian( &test_qm ); }
  void compare_gradient_and_hessian()
    { tester.compare_eval_with_grad_and_eval_with_hessian( &test_qm );  }
  void compare_analytical_and_numerical_hessians()
    { compare_analytical_and_numerical_hessians( &test_qm ); }
  void test_symmetric_hessian_diagonal()
    { tester.test_symmetric_Hessian_diagonal_blocks( &test_qm ); }
  void test_weighted_hessians()
    { compare_analytical_and_numerical_hessians( &weight_qm ); }
  void test_hessian_with_fixed_vertices()
    { tester.test_hessian_with_fixed_vertex( &center_qm ); }

  void compare_indices_and_diagonal()
    { tester.compare_eval_with_indices_and_eval_with_diagonal( &test_qm ); }
  void compare_gradient_and_diagonal()
    { tester.compare_eval_with_grad_and_eval_with_diagonal( &test_qm ); }
  void compare_analytical_and_numerical_diagonals()
    { compare_analytical_and_numerical_diagonals( &test_qm ); }
  void test_weighted_diagonals()
    { compare_analytical_and_numerical_diagonals( &weight_qm ); }
  void test_diagonal_with_fixed_vertices()
    { tester.test_diagonal_with_fixed_vertex( &center_qm ); }
  
  void test_inverse_mean_ratio_grad();
  void test_inverse_mean_ratio_hess();
  void test_inverse_mean_ratio_hess_diag();
  void regression_inverse_mean_ratio_grad();
  void regression_inverse_mean_ratio_hess();

    // Delcare specialized versions of the functions from
    // QualityMetricTester because we surface elements must
    // be handled differently.  For a surface element in the XY plane,
    // the finite difference approximations of the derivatives will
    // have non-zero values for derivatives wrt Z coordinates while the
    // analytical derivative calculations will return all derivatives
    // wrt Z coordiantes as zero.
  void get_nonideal_element( EntityTopology type, PatchData& pd );
  void compare_analytical_and_numerical_gradients( QualityMetric* qm );
  void compare_analytical_and_numerical_hessians( QualityMetric* qm );
  void compare_analytical_and_numerical_diagonals( QualityMetric* qm );
  void compare_analytical_and_numerical_gradients( QualityMetric* qm, PatchData&, int dim );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TMPQualityMetricTest, "TMPQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TMPQualityMetricTest, "Unit");

void TMPQualityMetricTest::test_get_evaluations()
{
  SamplePoints corners_and_edges( true, true, false, false );
  TMPQualityMetric edge_metric( &corners_and_edges, &ideal, &faux_2d_zero, &faux_3d_zero );
  
  tester.test_get_sample_evaluations( &zero_qm, &corners );
  tester.test_get_sample_evaluations( &edge_metric, &corners_and_edges );
}
  

void TMPQualityMetricTest::test_get_element_evaluations()
{
  SamplePoints edges( false, true, false, false );
  TMPQualityMetric edge_metric( &edges, &ideal, &faux_2d_zero, &faux_3d_zero );
  
  tester.test_get_in_element_evaluations( &zero_qm );
  tester.test_get_in_element_evaluations( &edge_metric );
}

static double col_dot_prod( MsqMatrix<2,2>& m )
  { return m(0,0) * m(0,1) + m(1,0) * m(1,1); }

void TMPQualityMetricTest::test_evaluate_2D()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &faux_2d_pi, &faux_3d_zero );
  
    // test with aligned elements
  faux_2d_pi.count = faux_3d_zero.count = 0;
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_pi.count );
  CPPUNIT_ASSERT_EQUAL( 0, faux_3d_zero.count );
  ASSERT_MATRICES_EQUAL( faux_2d_pi.last_W, faux_2d_pi.last_A, 1e-6 );
  
    // test that columns are orthogonal for ideal quad element
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, col_dot_prod(faux_2d_pi.last_A), 1e-6 );
  
    // test with an element rotated about X-axis
  faux_2d_pi.count = faux_3d_zero.count = 0;
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  // rotate by 90 degrees about X axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    Vector3D orig = pd.vertex_by_index(i);
    Vector3D newp( orig[0], -orig[2], orig[1] );
    pd.set_vertex_coordinates( newp, i, err );
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_pi.count );
  CPPUNIT_ASSERT_EQUAL( 0, faux_3d_zero.count );
  ASSERT_MATRICES_EQUAL( faux_2d_pi.last_W, faux_2d_pi.last_A, 1e-6 );
  
    // test with an element rotated about Y-axis
  faux_2d_pi.count = faux_3d_zero.count = 0;
  tester.get_ideal_element( TRIANGLE, true, pd );
  // rotate by -90 degrees about Y axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    Vector3D orig = pd.vertex_by_index(i);
    Vector3D newp( orig[2], orig[1], -orig[0] );
    pd.set_vertex_coordinates( newp, i, err );
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_pi.count );
  CPPUNIT_ASSERT_EQUAL( 0, faux_3d_zero.count );
  ASSERT_MATRICES_EQUAL( faux_2d_pi.last_W, faux_2d_pi.last_A, 1e-6 );
}  
 
  
void TMPQualityMetricTest::test_evaluate_3D()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &faux_2d_zero, &faux_3d_two );
  
    // test with aligned elements
  faux_2d_zero.count = faux_3d_two.count = 0;
  tester.get_ideal_element( HEXAHEDRON, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_3d_two.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 0, faux_2d_zero.count );
  CPPUNIT_ASSERT_EQUAL( 1, faux_3d_two.count );
  ASSERT_MATRICES_EQUAL( faux_3d_two.last_W, faux_3d_two.last_A, 1e-6 );
  
    // test that columns are orthogonal for ideal hex element
  MsqMatrix<3,3> A = faux_3d_two.last_A;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(1) % A.column(2), 1e-6 );
  
    // test with rotated element
  faux_2d_zero.count = faux_3d_two.count = 0;
  tester.get_ideal_element( TETRAHEDRON, true, pd );
  // rotate by 90-degrees about X axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    Vector3D orig = pd.vertex_by_index(i);
    Vector3D newp( orig[0], -orig[2], orig[1] );
    pd.set_vertex_coordinates( newp, i, err );
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_3d_two.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 0, faux_2d_zero.count );
  CPPUNIT_ASSERT_EQUAL( 1, faux_3d_two.count );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( det(faux_3d_two.last_W),
                                det(faux_3d_two.last_A),
                                1e-6 );
}


void TMPQualityMetricTest::test_evaluate_2D_weight()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &e_weight, &faux_2d_pi, &faux_3d_zero );
  
  tester.get_ideal_element( TRIANGLE, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value*e_weight.value, value, DBL_EPSILON );
}


void TMPQualityMetricTest::test_evaluate_3D_weight()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &e_weight, &faux_2d_zero, &faux_3d_two );
  
  tester.get_ideal_element( PRISM, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_3d_two.value*e_weight.value, value, DBL_EPSILON );
}

void TMPQualityMetricTest::test_2d_eval_ortho_quad( unsigned dim )
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  bool locs[4] = {false, false, false, false};
  locs[dim] = true;
  
  SamplePoints pts( locs[0], locs[1], locs[2], locs[3] );
  TMPQualityMetric m( &pts, &ideal, &faux_2d_zero, &faux_3d_zero );
  faux_2d_zero.count = faux_3d_zero.count = 0;
  
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_zero.count );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, col_dot_prod(faux_2d_zero.last_A), DBL_EPSILON );
}  

void TMPQualityMetricTest::test_3d_eval_ortho_hex( unsigned dim )
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  bool locs[4] = {false, false, false, false};
  locs[dim] = true;
  
  SamplePoints pts( locs[0], locs[1], locs[2], locs[3] );
  TMPQualityMetric m( &pts, &ideal, &faux_2d_zero, &faux_3d_zero );
  faux_2d_zero.count = faux_3d_zero.count = 0;
  
  tester.get_ideal_element( HEXAHEDRON, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( 1, faux_3d_zero.count );
  
    // test that columns are orthogonal for ideal hex element
  MsqMatrix<3,3> A = faux_3d_zero.last_A;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(1) % A.column(2), 1e-6 );
}  


void TMPQualityMetricTest::test_gradient_2D()
{
  MsqPrintError err(msq_stdio::cout);
  
    // check for expected value at center of flattened hex
  
    // construct flattened hex
  const double y = 0.5;
  const double vertices[] = { 0.0, 0.0, 0.0,
                              1.0, 0.0, 0.0,
                              1.0, y  , 0.0,
                              0.0, y  , 0.0 };
  size_t conn[8] = { 0, 1, 2, 3 };
  PatchData pd;
  pd.fill( 4, vertices, 1, QUADRILATERAL, conn, 0, err );
  ASSERT_NO_ERROR(err);
  
    // calculate Jacobian matrix at element center
  const double corner_xi[4][2] = { { -1, -1 },
                                   {  1, -1 },
                                   {  1,  1 },
                                   { -1,  1 } } ;
  MsqMatrix<4,2> coeff_derivs(&corner_xi[0][0]);
  coeff_derivs *= 0.25;  // derivatives of bilinear map at quad center
  MsqMatrix<4,3> coords( vertices );
  MsqMatrix<3,2> J = transpose(coords) * coeff_derivs;
    // calculate expected metric value
  const double expt_val = sqr_Frobenius( J );
    // calculate derivative for each element vertex
  MsqVector<3> expt_grad[4];
  for (int v = 0; v < 4; ++v)
    expt_grad[v] = 2 * J * transpose( coeff_derivs.row(v) );
    
  
    // construct metric
  pd.set_mapping_functions( &mf );
  SamplePoints center( false, false, true, false );
  TestGradTargetMetric2D tm;
  IdealTargetCalculator tc;
  TMPQualityMetric m( &center, &tc, &tm, 0 );
  PlanarDomain plane( PlanarDomain::XY );
  pd.set_domain( &plane );
  
    // evaluate metric
  double act_val;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> act_grad;
  m.evaluate_with_gradient( pd, 0, act_val, indices, act_grad, err );
  ASSERT_NO_ERROR(err);
  
    // compare values
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expt_val, act_val, 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[0]].data()), act_grad[0], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[1]].data()), act_grad[1], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[2]].data()), act_grad[2], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[3]].data()), act_grad[3], 1e-10 );

    // check numerical approx of gradient
  m.QualityMetric::evaluate_with_gradient( pd, 0, act_val, indices, act_grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expt_val, act_val, 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[0]].data()), act_grad[0], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[1]].data()), act_grad[1], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[2]].data()), act_grad[2], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[3]].data()), act_grad[3], 1e-5 );
}


void TMPQualityMetricTest::test_gradient_3D()
{
  MsqPrintError err(msq_stdio::cout);
  
    // check for expected value at center of flattened hex
  
    // construct flattened hex
  const double z = 0.5;
  const double vertices[] = { 0.0, 0.0, 0.0,
                              1.0, 0.0, 0.0,
                              1.0, 1.0, 0.0,
                              0.0, 1.0, 0.0,
                              0.0, 0.0, z,
                              1.0, 0.0, z,
                              1.0, 1.0, z,
                              0.0, 1.0, z };
  size_t conn[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  PatchData pd;
  pd.fill( 8, vertices, 1, HEXAHEDRON, conn, 0, err );
  ASSERT_NO_ERROR(err);
  
    // calculate Jacobian matrix at element center
  const double corner_xi[8][3] = { { -1, -1, -1 },
                                   {  1, -1, -1 },
                                   {  1,  1, -1 },
                                   { -1,  1, -1 },
                                   { -1, -1,  1 },
                                   {  1, -1,  1 },
                                   {  1,  1,  1 },
                                   { -1,  1,  1 } };
  MsqMatrix<8,3> coeff_derivs(&corner_xi[0][0]);
  coeff_derivs *= 0.125;  // derivatives of trilinear map at hex center
  MsqMatrix<8,3> coords( vertices );
  MsqMatrix<3,3> J = transpose(coords) * coeff_derivs;
    // calculate expected metric value
  const double expt_val = sqr_Frobenius( J );
    // calculate derivative for each element vertex
  MsqVector<3> expt_grad[8];
  for (int v = 0; v < 8; ++v)
    expt_grad[v] = 2 * J * transpose( coeff_derivs.row(v) );
    
  
    // construct metric
  pd.set_mapping_functions( &mf );
  SamplePoints center( false, false, false, true );
  TestGradTargetMetric3D tm;
  IdealTargetCalculator tc;
  TMPQualityMetric m( &center, &tc, 0, &tm );
  
    // evaluate metric
  double act_val;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> act_grad;
  m.evaluate_with_gradient( pd, 0, act_val, indices, act_grad, err );
  ASSERT_NO_ERROR(err);
  
    // compare values
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expt_val, act_val, 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[0]].data()), act_grad[0], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[1]].data()), act_grad[1], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[2]].data()), act_grad[2], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[3]].data()), act_grad[3], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[4]].data()), act_grad[4], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[5]].data()), act_grad[5], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[6]].data()), act_grad[6], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[7]].data()), act_grad[7], 1e-10 );

    // check numerical approx of gradient
  m.QualityMetric::evaluate_with_gradient( pd, 0, act_val, indices, act_grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expt_val, act_val, 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[0]].data()), act_grad[0], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[1]].data()), act_grad[1], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[2]].data()), act_grad[2], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[3]].data()), act_grad[3], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[4]].data()), act_grad[4], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[5]].data()), act_grad[5], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[6]].data()), act_grad[6], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[7]].data()), act_grad[7], 1e-5 );
}

void TMPQualityMetricTest::get_nonideal_element( EntityTopology type,
                                                 PatchData& pd )
{
  tester.get_nonideal_element( type, pd, true );
    // Callers assume surface elements are in XY plane.
    // Verify this assumption.
  if (TopologyInfo::dimension(type) == 2) {
    for (size_t i = 0; i < pd.num_nodes(); ++i) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL( pd.vertex_by_index(i)[2], 0.0, 1e-6 );
    }
  }
}

    // Delcare specialized versions of the functions from
    // QualityMetricTester because we surface elements must
    // be handled differently.  For a surface element in the XY plane,
    // the finite difference approximations of the derivatives will
    // have non-zero values for derivatives wrt Z coordinates while the
    // analytical derivative calculations will return all derivatives
    // wrt Z coordiantes as zero.
void TMPQualityMetricTest::compare_analytical_and_numerical_gradients( QualityMetric* qm )
{
  PatchData pd;
  const EntityTopology types[] = { TRIANGLE,
                                   QUADRILATERAL,
                                   TETRAHEDRON,
                                   PYRAMID,
                                   PRISM,
                                   HEXAHEDRON };
  const int num_types = sizeof(types)/sizeof(types[0]);
  for (int i = 0; i < num_types; ++i) {
    get_nonideal_element( types[i], pd );
    compare_analytical_and_numerical_gradients( qm, pd, TopologyInfo::dimension(types[i]) );
  }
}

void TMPQualityMetricTest::compare_analytical_and_numerical_gradients( 
                                                      QualityMetric* qm,
                                                      PatchData& pd,
                                                      int dim )
{
  MsqPrintError err( msq_stdio::cout );

  msq_std::vector<size_t> handles, indices1, indices2;
  msq_std::vector<Vector3D> grad1, grad2;
  double qm_val1, qm_val2;
  bool rval;

  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->QualityMetric::evaluate_with_gradient( pd, handles[j], qm_val1, indices1, grad1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

      // For analytical gradient of a 2D element in the XY plane, 
      // we expect all Z terms to be zero.
    if (dim == 2)
      for (size_t k = 0; k < grad1.size(); ++k)
        grad1[k][2] = 0.0; 

    rval = qm->evaluate_with_gradient( pd, handles[j], qm_val2, indices2, grad2, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );
    CPPUNIT_ASSERT( !indices1.empty() );

    msq_std::vector<size_t>::iterator it1, it2;
    for (it1 = indices1.begin(); it1 != indices1.end(); ++it1) {
      it2 = msq_std::find( indices2.begin(), indices2.end(), *it1 );
      CPPUNIT_ASSERT( it2 != indices2.end() );

      size_t idx1 = it1 - indices1.begin();
      size_t idx2 = it2 - indices2.begin();
      CPPUNIT_ASSERT_VECTORS_EQUAL( grad1[idx1], grad2[idx2], 0.01 );
    }
  }
}


    // Delcare specialized versions of the functions from
    // QualityMetricTester because we surface elements must
    // be handled differently.  For a surface element in the XY plane,
    // the finite difference approximations of the derivatives will
    // have non-zero values for derivatives wrt Z coordinates while the
    // analytical derivative calculations will return all derivatives
    // wrt Z coordiantes as zero.
void TMPQualityMetricTest::compare_analytical_and_numerical_hessians( QualityMetric* qm )
{
  MsqPrintError err( msq_stdio::cout );
  PatchData pd;
  const EntityTopology types[] = { TRIANGLE,
                                   QUADRILATERAL,
                                   TETRAHEDRON,
                                   PYRAMID,
                                   PRISM,
                                   HEXAHEDRON };
  const int num_types = sizeof(types)/sizeof(types[0]);
  for (int i = 0; i < num_types; ++i) {
    get_nonideal_element( types[i], pd );

    msq_std::vector<size_t> handles, indices1, indices2;
    msq_std::vector<Vector3D> grad1, grad2;
    msq_std::vector<Matrix3D> Hess1, Hess2;
    double qm_val1, qm_val2;
    bool rval;

    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    for (size_t j = 0; j < handles.size(); ++j) {
      rval = qm->QualityMetric::evaluate_with_Hessian( pd, handles[j], qm_val1, indices1, grad1, Hess1, err );
      CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      CPPUNIT_ASSERT( rval );
       
        // For analytical gradient of a 2D element in the XY plane, 
        // we expect all Z terms to be zero.
      if (TopologyInfo::dimension(types[i]) == 2)
        for (size_t k = 0; k < Hess1.size(); ++k) 
          Hess1[k](0,2) = Hess1[k](1,2) = Hess1[k](2,0) 
            = Hess1[k](2,1) = Hess1[k](2,2) = 0.0;

      rval = qm->evaluate_with_Hessian( pd, handles[j], qm_val2, indices2, grad2, Hess2, err );
      CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      CPPUNIT_ASSERT( rval );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
      CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );
      CPPUNIT_ASSERT( !indices1.empty() );

      msq_std::vector<size_t>::iterator it;
      unsigned h = 0;
      for (unsigned r = 0; r < indices1.size(); ++r) {
        it = msq_std::find( indices2.begin(), indices2.end(), indices1[r] );
        CPPUNIT_ASSERT( it != indices2.end() );
        unsigned r2 = it - indices2.begin();

        for (unsigned c = r; c < indices1.size(); ++c, ++h) {
          it = msq_std::find( indices2.begin(), indices2.end(), indices1[c] );
          CPPUNIT_ASSERT( it != indices2.end() );
          unsigned c2 = it - indices2.begin();

          unsigned h2;
          if (r2 <= c2) 
            h2 = indices2.size()*r - r*(r+1)/2 + c;
          else
            h2 = indices2.size()*c - c*(c+1)/2 + r;

          //if (!utest_mat_equal(Hess1[h],Hess2[h2],0.001))
          //  assert(false);
          CPPUNIT_ASSERT_MATRICES_EQUAL( Hess1[h], Hess2[h2], 0.001 );
        }
      }
    }
  }
}

    // Delcare specialized versions of the functions from
    // QualityMetricTester because we surface elements must
    // be handled differently.  For a surface element in the XY plane,
    // the finite difference approximations of the derivatives will
    // have non-zero values for derivatives wrt Z coordinates while the
    // analytical derivative calculations will return all derivatives
    // wrt Z coordiantes as zero.
void TMPQualityMetricTest::compare_analytical_and_numerical_diagonals( QualityMetric* qm )
{
  MsqPrintError err( msq_stdio::cout );
  PatchData pd;
  const EntityTopology types[] = { TRIANGLE,
                                   QUADRILATERAL,
                                   TETRAHEDRON,
                                   PYRAMID,
                                   PRISM,
                                   HEXAHEDRON };
  const int num_types = sizeof(types)/sizeof(types[0]);
  for (int i = 0; i < num_types; ++i) {
    get_nonideal_element( types[i], pd );

    msq_std::vector<size_t> handles, indices1, indices2;
    msq_std::vector<Vector3D> grad1, grad2;
    msq_std::vector<Matrix3D> Hess1;
    msq_std::vector<SymMatrix3D> Hess2;
    double qm_val1, qm_val2;
    bool rval;

    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    for (size_t j = 0; j < handles.size(); ++j) {
      rval = qm->QualityMetric::evaluate_with_Hessian( pd, handles[j], qm_val1, indices1, grad1, Hess1, err );
      CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      CPPUNIT_ASSERT( rval );
       
        // For analytical gradient of a 2D element in the XY plane, 
        // we expect all Z terms to be zero.
      if (TopologyInfo::dimension(types[i]) == 2)
        for (size_t k = 0; k < Hess1.size(); ++k) 
          Hess1[k](0,2) = Hess1[k](1,2) = Hess1[k](2,0) 
            = Hess1[k](2,1) = Hess1[k](2,2) = 0.0;

      rval = qm->evaluate_with_Hessian_diagonal( pd, handles[j], qm_val2, indices2, grad2, Hess2, err );
      CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      CPPUNIT_ASSERT( rval );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
      CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );
      CPPUNIT_ASSERT( !indices1.empty() );
      CPPUNIT_ASSERT_EQUAL( indices1.size() * (indices1.size()+1) / 2, Hess1.size() );
      CPPUNIT_ASSERT_EQUAL( indices2.size(), Hess2.size() );

      size_t h = 0;
      msq_std::vector<size_t>::iterator it;
      for (unsigned r = 0; r < indices1.size(); ++r) {
        it = msq_std::find( indices2.begin(), indices2.end(), indices1[r] );
        CPPUNIT_ASSERT( it != indices2.end() );
        unsigned r2 = it - indices2.begin();
        //if (!utest_mat_equal(Hess1[h],Hess2[r2],0.001))
        //  assert(false);
        CPPUNIT_ASSERT_MATRICES_EQUAL( Hess1[h], Hess2[r2], 0.001 );
        h += indices1.size() - r;
      }
    }
  }
}

  
void TMPQualityMetricTest::test_inverse_mean_ratio_grad()
{
  InverseMeanRatio2D tm_2d;
  InverseMeanRatio3D tm_3d;
  SamplePoints pts(true);
  IdealTargetCalculator target(false);
  TMPQualityMetric metric( &pts, &target, &tm_2d, &tm_3d );
  ElementPMeanP avg( 1.0, &metric );
  
  tester.test_gradient_reflects_quality( &metric );
  compare_analytical_and_numerical_gradients( &metric );
  tester.test_gradient_with_fixed_vertex( &avg );
}


void TMPQualityMetricTest::test_inverse_mean_ratio_hess()
{
  InverseMeanRatio2D tm_2d;
  InverseMeanRatio3D tm_3d;
  SamplePoints pts(true);
  IdealTargetCalculator target(false);
  TMPQualityMetric metric( &pts, &target, &tm_2d, &tm_3d );
  ElementPMeanP avg( 1.0, &metric );
 
  compare_analytical_and_numerical_hessians( &metric );
  tester.test_symmetric_Hessian_diagonal_blocks( &metric );
  tester.test_hessian_with_fixed_vertex( &avg );
}

void TMPQualityMetricTest::test_inverse_mean_ratio_hess_diag()
{
  InverseMeanRatio2D tm_2d;
  InverseMeanRatio3D tm_3d;
  SamplePoints pts(true);
  IdealTargetCalculator target(false);
  TMPQualityMetric metric( &pts, &target, &tm_2d, &tm_3d );
  
  compare_analytical_and_numerical_diagonals( &metric );
  tester.compare_eval_with_diag_and_eval_with_hessian( &metric );
}
  
void TMPQualityMetricTest::regression_inverse_mean_ratio_grad()
{
  MsqError err;
  InverseMeanRatio2D tm_2d;
  SamplePoints pts(true);
  IdealTargetCalculator target(false);
  TMPQualityMetric metric( &pts, &target, &tm_2d, 0 );
  const double coords[] = { -0.80000000000000004, -0.80000000000000004, 0,
                             0.00000000000000000,  2.00000000000000000, 0,
                            -1.73205079999999990,  1.00000000000000000, 0 };
  const size_t indices[] = { 0, 1, 2 };
  PatchData pd;
  pd.fill( 3, coords, 1, TRIANGLE, indices, 0, err );
  pd.set_mapping_functions( &mf );
  PlanarDomain dom( PlanarDomain::XY, coords[0] );
  pd.set_domain( &dom );
  
  IdealWeightInverseMeanRatio ref_metric;
  
  double exp_val, act_val;
  std::vector<size_t> exp_idx, act_idx, handles;
  std::vector<Vector3D> exp_grad, act_grad;
  
  handles.clear();
  ref_metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  const size_t hand = handles.front();
  handles.clear();
  metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  CPPUNIT_ASSERT_EQUAL( hand, handles.front() );
  
  bool exp_rval, act_rval;
  exp_rval = ref_metric.evaluate_with_gradient( pd, hand, exp_val, exp_idx, exp_grad, err );
  ASSERT_NO_ERROR(err);
  act_rval = metric.evaluate_with_gradient( pd, hand, act_val, act_idx, act_grad, err );
  ASSERT_NO_ERROR(err);
  
  CPPUNIT_ASSERT( exp_rval );
  CPPUNIT_ASSERT( act_rval );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_val, act_val, 1e-5 );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, exp_idx.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, act_idx.size() );
  
  std::vector<size_t> sorted(exp_idx);
  msq_std::sort( sorted.begin(), sorted.end() );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, sorted[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, sorted[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, sorted[2] );
  
  sorted = act_idx;
  msq_std::sort( sorted.begin(), sorted.end() );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, sorted[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, sorted[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, sorted[2] );
  
  const size_t idx_map[] = { 
    msq_std::find(act_idx.begin(),act_idx.end(),exp_idx[0]) - act_idx.begin(),
    msq_std::find(act_idx.begin(),act_idx.end(),exp_idx[1]) - act_idx.begin(),
    msq_std::find(act_idx.begin(),act_idx.end(),exp_idx[2]) - act_idx.begin() };
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[idx_map[0]], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[idx_map[1]], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[2], act_grad[idx_map[2]], 1e-5 );
}

  
void TMPQualityMetricTest::regression_inverse_mean_ratio_hess()
{
  MsqError err;
  InverseMeanRatio2D tm_2d;
  SamplePoints pts(true);
  IdealTargetCalculator target(false);
  TMPQualityMetric metric( &pts, &target, &tm_2d, 0 );
  const double coords[] = { 4.158984727, 4.6570859130000004, 5,
                            4.51742825, 4.51742825, 5,
                            4.3103448279999999, 5, 5 };
  const bool fixed[] = { false, false, true };
  const size_t indices[] = { 0, 1, 2 };
  PatchData pd;
  pd.fill( 3, coords, 1, TRIANGLE, indices, fixed, err );
  pd.set_mapping_functions( &mf );
  PlanarDomain dom( PlanarDomain::XY, coords[2] );
  pd.set_domain( &dom );
  
  IdealWeightInverseMeanRatio ref_metric;
  
  double exp_val, act_val;
  std::vector<size_t> exp_idx, act_idx, handles;
  std::vector<Vector3D> exp_grad, act_grad;
  std::vector<Matrix3D> exp_hess, act_hess;
  
  handles.clear();
  ref_metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  const size_t hand = handles.front();
  handles.clear();
  metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  CPPUNIT_ASSERT_EQUAL( hand, handles.front() );
  
  bool exp_rval, act_rval;
  exp_rval = ref_metric.evaluate_with_Hessian( pd, hand, exp_val, exp_idx, exp_grad, exp_hess, err );
  ASSERT_NO_ERROR(err);
  act_rval = metric.evaluate_with_Hessian( pd, hand, act_val, act_idx, act_grad, act_hess, err );
  ASSERT_NO_ERROR(err);
  
  CPPUNIT_ASSERT( exp_rval );
  CPPUNIT_ASSERT( act_rval );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_val, act_val, 1e-5 );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, exp_idx.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, act_idx.size() );

  // zero derivatives with respect to Z
  for (int i = 0; i < 3; ++i) {
    exp_grad[i][2] = 0.0;
    for (int j = 0; j < 3; ++j) 
      exp_hess[i][j][2] = exp_hess[i][2][j] = 0.0;
  }
  
  if (act_idx[0] == exp_idx[0]) {
    CPPUNIT_ASSERT_EQUAL( exp_idx[1], act_idx[1] );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[0], 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[1], 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], act_hess[0], 1e-4 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[1], act_hess[1], 1e-4 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[2], act_hess[2], 1e-4 );
  }
  else {
    CPPUNIT_ASSERT_EQUAL( exp_idx[0], act_idx[1] );
    CPPUNIT_ASSERT_EQUAL( exp_idx[1], act_idx[0] );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[1], 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[0], 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], act_hess[2], 1e-4 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[1], transpose(act_hess[1]), 1e-4 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[2], act_hess[0], 1e-4 );
  }
}
