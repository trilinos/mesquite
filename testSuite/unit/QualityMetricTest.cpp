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
   
  ***************************************************************** */

/*! \file QualityMetricTest.cpp

Unit testing for the QualityMetric class
\author Jasno Kraftcheck
*/
#include "Mesquite.hpp"
#include "VertexQM.hpp"
#include "ElementQM.hpp"
#include "IdealElements.hpp"
#include "UnitUtil.hpp"
#include "TopologyInfo.hpp"
#include "PatchData.hpp"
#include "cppunit/extensions/HelperMacros.h"

#include <string>

using namespace Mesquite;

class QualityMetricTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QualityMetricTest);
  CPPUNIT_TEST (test_fixed_vertex_list);
  CPPUNIT_TEST (test_remove_fixed_gradients);
  CPPUNIT_TEST (test_remove_fixed_hessians);
  CPPUNIT_TEST( test_gradient_constant );
  CPPUNIT_TEST( test_gradient_linear );
  CPPUNIT_TEST( test_gradient_parabolic );
  CPPUNIT_TEST( test_Hessian_constant );
  CPPUNIT_TEST( test_Hessian_linear );
  CPPUNIT_TEST( test_Hessian_parabolic );
  CPPUNIT_TEST_SUITE_END();
  PatchData tri_pd;

public:
  void setUp();
  
  void test_fixed_vertex_list();
  void test_remove_fixed_gradients();
  void test_remove_fixed_hessians();
  void test_gradient_constant();
  void test_gradient_linear();
  void test_gradient_parabolic();
  void test_Hessian_constant();
  void test_Hessian_linear();
  void test_Hessian_parabolic();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "QualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "Unit");


// Create a single-triangle patch with one fixed vertex
void QualityMetricTest::setUp()
{
  MsqError err;
  const double vtx_coords[] = { 1, 0, 0,
                                0, 1, 0,
                                0, 0, 0 };
  const size_t connectivity[] = { 0, 1, 2 };
  const bool fixed[] = { false, false, true };
  tri_pd.fill( 3, vtx_coords, 1, TRIANGLE, connectivity, fixed, err );
  CPPUNIT_ASSERT(!err);
}

// tolerance on numerical gradient/hessian values
const double EPSILON = 1e-4;

/**\brief Fake quality metric for testing numerial gradient
 *
 * Implement a vertex metric for which the "quality" at a given
 * vertex is the value of the x-coordinate of that vertex.  Thus
 * the gradient of the "quality" at every vertex should be {1,0,0}.
 */
class LinearVertexMetric : public VertexQM
{
public:
  msq_std::string get_name() const { return "Fake metric for testing numerical gradient"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t vtx_idx, double& value, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[0];
    return true;
  }
  bool evaluate_with_indices( PatchData& pd, size_t vtx_idx, double& value, msq_std::vector<size_t>& indices, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[0];
    indices.resize(1);
    indices[0] = vtx_idx;
    return true;
  }
};

/**\brief Fake quality metric for testing numerial gradient
 *
 * Implement an element metric for which the "quality" is always 1.0.
 * The resulting gradient and Hessian should always be zero.
 */
class ConstantElementMetric : public ElementQM
{
public:
  msq_std::string get_name() const { return "Fake metric for testing numerical gradient"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t , double& value, MsqError& )
    { value = 1.0; return true; }
  bool evaluate_with_indices( PatchData& pd, size_t elem_idx, double& value, msq_std::vector<size_t>& indices, MsqError& )
  { 
    MsqMeshEntity& elem = pd.element_by_index( elem_idx );
    unsigned nv = elem.node_count();
    const size_t* conn = elem.get_vertex_index_array();
    indices.clear();
    for (unsigned i = 0; i < nv; ++i)
      if (conn[i] < pd.num_free_vertices())
        indices.push_back( conn[i] );

    value = 1.0; 
    return true; 
  }
};

/**\brief Fake quality metric for testing numerial gradient
 *
 * Implement a vertex metric for which the "quality" at a given
 * vertex is the square of the value of the y-coordinate of that vertex.  Thus
 * the Hessian of the "quality" at every vertex should be {2,0,0}.
 */
class ParabolicVertexMetric : public VertexQM
{
public:
  msq_std::string get_name() const { return "Fake metric for testing numerical gradient"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t vtx_idx, double& value, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[1];
    value *= value;
    return true;
  }
  bool evaluate_with_indices( PatchData& pd, size_t vtx_idx, double& value, msq_std::vector<size_t>& indices, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[1];
    value *= value;
    indices.resize(1);
    indices[0] = vtx_idx;
    return true;
  }
};

void QualityMetricTest::test_fixed_vertex_list()
{
  // define a patch of four quads such that
  // the number of fixed vertices in each quad
  // varies from 0 to 3
  /*   2------1-----8
   *   |      |     |\
   *   | (0)  | (3) | \
   *   |      |     |  \
   *   3------0-----7   \
   *   |      |     |\   \
   *   | (1)  | (2) | \   \
   *   |      |     |  \   \
   *   4------5-----6   \   \
   *           \____\____\___\__ fixed
   */
  const double coords[] = { 0, 0, 0,
                            0, 1, 0,
                           -1, 1, 0,
                           -1, 0, 0,
                           -1,-1, 0,
                            0,-1, 0,
                            1,-1, 0,
                            1, 0, 0,
                            1, 1, 0 };
  const size_t conn[] = { 0, 1, 2, 3,
                          3, 4, 5, 0, 
                          6, 7, 0, 5,
                          0, 7, 8, 1 };
  const bool fixed[] = { false, false, false, 
                         false, false, true,
                         true,  true,  true };
  
  MsqPrintError err(msq_stdio::cout);
  PatchData pd;
  pd.fill( 9, coords, 4, QUADRILATERAL, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  uint32_t bits;
  msq_std::vector<size_t> indices;
  msq_std::vector<size_t>::iterator it;
  const size_t* verts;
  unsigned i;
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(0), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)4, indices.size() );
  CPPUNIT_ASSERT_EQUAL( (uint32_t)0, bits&0xF );
  CPPUNIT_ASSERT( pd.num_free_vertices() > *msq_std::max_element(indices.begin(), indices.end()) );
  verts = pd.element_by_index(0).get_vertex_index_array();
  for (i = 0; i < 4; ++i) 
    CPPUNIT_ASSERT( msq_stdc::find( verts, verts+4, indices[i] ) != verts+4 );
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(1), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, indices.size() );
  verts = pd.element_by_index(1).get_vertex_index_array();
  for (i = 0; i < 4; ++i) {
    it = msq_std::find( indices.begin(), indices.end(), verts[i] );
    if (verts[i] < pd.num_free_vertices()) {
      CPPUNIT_ASSERT( it != indices.end() );
      CPPUNIT_ASSERT_EQUAL( 0u, bits & (1<<i) );
    }
    else {
      CPPUNIT_ASSERT( it == indices.end() );
      CPPUNIT_ASSERT_EQUAL( 1u, (bits>>i) & 1 );
    }
  }
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(2), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() );
  verts = pd.element_by_index(2).get_vertex_index_array();
  for (i = 0; i < 4; ++i) {
    it = msq_std::find( indices.begin(), indices.end(), verts[i] );
    if (verts[i] < pd.num_free_vertices()) {
      CPPUNIT_ASSERT( it != indices.end() );
      CPPUNIT_ASSERT_EQUAL( 0u, bits & (1<<i) );
    }
    else {
      CPPUNIT_ASSERT( it == indices.end() );
      CPPUNIT_ASSERT_EQUAL( 1u, (bits>>i) & 1 );
    }
  }
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(3), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, indices.size() );
  verts = pd.element_by_index(3).get_vertex_index_array();
  for (i = 0; i < 4; ++i) {
    it = msq_std::find( indices.begin(), indices.end(), verts[i] );
    if (verts[i] < pd.num_free_vertices()) {
      CPPUNIT_ASSERT( it != indices.end() );
      CPPUNIT_ASSERT_EQUAL( 0u, bits & (1<<i) );
    }
    else {
      CPPUNIT_ASSERT( it == indices.end() );
      CPPUNIT_ASSERT_EQUAL( 1u, (bits>>i) & 1 );
    }
  }
}  

void QualityMetricTest::test_remove_fixed_gradients()
{
    // define a list of vectors
  msq_std::vector<Vector3D> grads(6);
  grads[0] = Vector3D(0,0,0);
  grads[1] = Vector3D(1,1,1);
  grads[2] = Vector3D(2,2,2);
  grads[3] = Vector3D(3,3,3);
  grads[4] = Vector3D(4,4,4);
  grads[5] = Vector3D(5,5,5);
    // remove the first, third, and last
  uint32_t flags = 1u | 4u | 32u;
    // call function, choose element type w/ correct number of corners
  QualityMetric::remove_fixed_gradients( PRISM, flags, grads );
    // check results
  CPPUNIT_ASSERT_EQUAL( (size_t)3, grads.size() );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1,1,1), grads[0], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(3,3,3), grads[1], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(4,4,4), grads[2], DBL_EPSILON );
}

void QualityMetricTest::test_remove_fixed_hessians()
{
    // define Hessian matrix for a quadrilateral
  Matrix3D m[10] = {
    Matrix3D(0.0), Matrix3D(1.0), Matrix3D(2.0), Matrix3D(3.0),
                   Matrix3D(4.0), Matrix3D(5.0), Matrix3D(6.0),
                                  Matrix3D(7.0), Matrix3D(8.0),
                                                 Matrix3D(9.0)
  };
    // convert to std::vector
  msq_std::vector<Matrix3D> hess(10);
  msq_std::copy( m, m+10, hess.begin() );
    // mark fist and third vertices as fixed
  uint32_t flags = 1u | 4u;
    // call function to remove grads for fixed vertices
  QualityMetric::remove_fixed_hessians( QUADRILATERAL, flags, hess );
    // the submatrix with the first and third rows/columns removed should be
    // { 4, 6,
    //      9 }
  CPPUNIT_ASSERT_EQUAL( (size_t)3, hess.size() );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(4.0), hess[0], DBL_EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(6.0), hess[1], DBL_EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(9.0), hess[2], DBL_EPSILON );
}


void QualityMetricTest::test_gradient_constant()
{
  MsqError err;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> gradient;
  double value;
  size_t ELEMENT = 0;

  ConstantElementMetric constant;
  QualityMetric& qm = constant;

  qm.evaluate_with_gradient( tri_pd, ELEMENT, value, indices, gradient, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(indices.size(), gradient.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());  // two free vertices in triangle
  
  const Vector3D& g1 = gradient[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, g1[0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, g1[1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, g1[2], EPSILON );
  
  const Vector3D& g2 = gradient[1];
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, g2[0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, g2[1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, g2[2], EPSILON );
}


void QualityMetricTest::test_gradient_linear()
{
  MsqError err;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> gradient;
  double value;
  size_t VERTEX = 0;
  
  LinearVertexMetric linear;
  QualityMetric& qm = linear;

  qm.evaluate_with_gradient( tri_pd, VERTEX, value, indices, gradient, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(indices.size(), gradient.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  
  const Vector3D& grad = gradient[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, grad[0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, grad[1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, grad[2], EPSILON );
}


void QualityMetricTest::test_gradient_parabolic()
{
  const size_t VERTEX = 1;  // pick vertex with non-zero Y-coordinate
  MsqError err;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> gradient;
  double value;

  ParabolicVertexMetric parab;
  QualityMetric& qm = parab;

  qm.evaluate_with_gradient( tri_pd, VERTEX, value, indices, gradient, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(indices.size(), gradient.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());  // two free vertices in triangle
  
  const double expected_y = 2*tri_pd.vertex_by_index(VERTEX)[1];
  const Vector3D& g = gradient[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL(        0.0, g[0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expected_y, g[1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(        0.0, g[2], EPSILON );
}


void QualityMetricTest::test_Hessian_constant()
{
  MsqError err;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> gradient;
  msq_std::vector<Matrix3D> Hessian;
  double value;
  size_t ELEMENT = 0;

  ConstantElementMetric constant;
  QualityMetric& qm = constant;
  
  qm.evaluate_with_Hessian( tri_pd, ELEMENT, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(indices.size()*(indices.size()+1)/2, Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());
}


void QualityMetricTest::test_Hessian_linear()
{
  MsqError err;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> gradient;
  msq_std::vector<Matrix3D> Hessian;
  double value;
  size_t VERTEX = 0;
  
  LinearVertexMetric linear;
  QualityMetric& qm = linear;

  qm.evaluate_with_Hessian( tri_pd, VERTEX, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(indices.size()*(indices.size()+1)/2, Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  
  const Matrix3D& H = Hessian[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[0][0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[0][1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[0][2], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[1][0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[1][1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[1][2], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[2][0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[2][1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[2][2], EPSILON );
}


void QualityMetricTest::test_Hessian_parabolic()
{
  MsqError err;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> gradient;
  msq_std::vector<Matrix3D> Hessian;
  double value;
  size_t VERTEX = 0;

  ParabolicVertexMetric parab;
  QualityMetric& qm = parab;
  
  qm.evaluate_with_Hessian( tri_pd, VERTEX, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(indices.size()*(indices.size()+1)/2, Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  
  const Matrix3D& H = Hessian[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[0][0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[0][1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[0][2], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[1][0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, H[1][1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[1][2], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[2][0], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[2][1], EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, H[2][2], EPSILON );
}

