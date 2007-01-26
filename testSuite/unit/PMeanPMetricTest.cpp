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


#include "ElementPMeanP.hpp"
#include "VertexPMeanP.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "PatchDataInstances.hpp"
#include "UnitUtil.hpp"
#include "ElemSampleQM.hpp"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class PMeanPMetricTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(PMeanPMetricTest);
  CPPUNIT_TEST(test_get_metric_type);
  CPPUNIT_TEST(test_get_element_evaluations);
  CPPUNIT_TEST(test_get_vertex_evaluations);
  CPPUNIT_TEST(test_element_evaluate);
  CPPUNIT_TEST(test_vertex_evaluate);
  CPPUNIT_TEST(test_indices);
  CPPUNIT_TEST(test_gradient);
  CPPUNIT_TEST_SUITE_END();

  PatchData pd;

public:
  void setUp();

  void test_get_metric_type();
  void test_get_element_evaluations();
  void test_get_vertex_evaluations();
  void test_element_evaluate();
  void test_vertex_evaluate();
  void test_indices();
  void test_gradient();
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PMeanPMetricTest, "PMeanPMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PMeanPMetricTest, "Unit");

void PMeanPMetricTest::setUp() {
  MsqError err;
  create_four_quads_patch( pd, err );
  CPPUNIT_ASSERT(!err);
}

class FauxMetric : public ElemSampleQM
{
public:
  MetricType get_metric_type() const { return ELEMENT_BASED; }
  msq_std::string get_name() const { return "FauxMetrix"; }
  int get_negate_flag() const { return 1; }
  void get_evaluations( PatchData& pd, msq_std::vector<size_t>& h, bool free, MsqError& );
  void get_element_evaluations( PatchData&, size_t, msq_std::vector<size_t>&, MsqError& );
  bool evaluate( PatchData& pd, size_t h, double& v, MsqError& err );
  bool evaluate_with_indices( PatchData& pd, size_t h, double& v, 
                              msq_std::vector<size_t>& indices, MsqError& err );
  bool evaluate_with_gradient( PatchData& pd, size_t h, double& v, 
                              msq_std::vector<size_t>& indices, 
                              msq_std::vector<Vector3D>& grads,
                              MsqError& err );
};

void FauxMetric::get_evaluations( PatchData& pd, msq_std::vector<size_t>& h, bool free, MsqError& err )
{
  h.clear();
  for (size_t i = 0; i < pd.num_elements(); ++i)
    get_element_evaluations( pd, i, h, err );
}

void FauxMetric::get_element_evaluations( PatchData& pd, size_t h, msq_std::vector<size_t>& list, MsqError& )
{
  MsqMeshEntity& elem = pd.element_by_index(h);
  for (unsigned i = 0; i  < elem.corner_count(); ++i)
    list.push_back( handle( i, h ) );
}

bool FauxMetric::evaluate( PatchData& pd, size_t h, double& v, MsqError&  )
{
  size_t e = ElemSampleQM::elem( h );
  unsigned s = ElemSampleQM::sample( h );
  size_t* verts = pd.element_by_index(e).get_vertex_index_array();
  v = (double)(verts[s]);
  return true;
}

bool FauxMetric::evaluate_with_indices( PatchData& pd, size_t h, double& v, 
                              msq_std::vector<size_t>& indices, MsqError& err )
{
  evaluate( pd, h, v, err );
  indices.resize(3);
  size_t e = ElemSampleQM::elem( h );
  unsigned s = ElemSampleQM::sample( h );
  size_t* verts = pd.element_by_index(e).get_vertex_index_array();
  size_t n = pd.element_by_index(e).vertex_count();
  indices[0] = verts[s];
  indices[1] = verts[(s+1)%n];
  indices[2] = verts[(s+n-1)%n];
  return true;
}

bool FauxMetric::evaluate_with_gradient( PatchData& pd, size_t h, double& v, 
                              msq_std::vector<size_t>& indices, 
                              msq_std::vector<Vector3D>& grads,
                              MsqError& err )
{
  evaluate_with_indices( pd, h, v, indices, err );
  grads.clear();
  for (unsigned i = 0; i < indices.size(); ++i)
    grads.push_back( Vector3D( (double)indices[i], 0, 1 ) );
  return true;
}


void PMeanPMetricTest::test_get_metric_type()
{
  FauxMetric m;
  ElementPMeanP e( 1.0, &m );
  VertexPMeanP v( 1.0, &m );
  CPPUNIT_ASSERT( QualityMetric::ELEMENT_BASED == e.get_metric_type() );
  CPPUNIT_ASSERT( QualityMetric::VERTEX_BASED == v.get_metric_type() );
}

void PMeanPMetricTest::test_get_element_evaluations()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP e( 1.0, &m );
  msq_std::vector<size_t> handles;
    // test that handles array contains all elements
  e.get_evaluations( pd, handles, false, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
    // test that handles array contains all elements
  e.get_evaluations( pd, handles, true, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
}
  

void PMeanPMetricTest::test_get_vertex_evaluations()
{
  MsqError err;
  FauxMetric m;
  VertexPMeanP e( 1.0, &m );
  msq_std::vector<size_t> handles;
    // test that handles array contains all vertices
  e.get_evaluations( pd, handles, false, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_nodes(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
    // test that handles array contains all vertices
  e.get_evaluations( pd, handles, true, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_nodes(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
}

void PMeanPMetricTest::test_vertex_evaluate()
{
  MsqError err;
  FauxMetric m;
  VertexPMeanP m1( 1.0, &m );
  VertexPMeanP m2( 2.0, &m );
  
    // evaluate average around vertex
  double v1, v2;
  m1.evaluate( pd, 0, v1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate( pd, 0, v2, err );
  CPPUNIT_ASSERT(!err);

    // get elements adjacent to vertex
  size_t num_elem;
  const size_t* elems = pd.get_vertex_element_adjacencies( 0, num_elem, err );
  CPPUNIT_ASSERT(!err);
  
    // calculate expected values from underyling metric
  double ev1 = 0.0, ev2 = 0.0;
  for (unsigned i = 0; i < num_elem; ++i) {
    const MsqMeshEntity& elem = pd.element_by_index( elems[i] );
    const size_t* verts = elem.get_vertex_index_array();
    const size_t* end = verts + elem.node_count();
    const size_t* p = msq_std::find( verts, end, (size_t)0 );
    CPPUNIT_ASSERT( p < end );
    size_t h = ElemSampleQM::handle( p - verts, elems[i] );
  
    double v;
    m.evaluate( pd, h, v, err );
    CPPUNIT_ASSERT(!err);
    ev1 += v;
    ev2 += v*v;
  }
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev1, v1, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev2, v2, 1e-6 );
}

void PMeanPMetricTest::test_element_evaluate()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m1( 1.0, &m );
  ElementPMeanP m2( 2.0, &m );
  
    // evaluate average over element
  double v1, v2;
  m1.evaluate( pd, 0, v1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate( pd, 0, v2, err );
  CPPUNIT_ASSERT(!err);
  
    // get sample points within element
  msq_std::vector<size_t> handles;
  m.get_element_evaluations( pd, 0, handles, err );
  CPPUNIT_ASSERT(!err);
  
    // calculate expected values from underyling metric
  double ev1 = 0.0, ev2 = 0.0;
  for (unsigned i = 0; i < handles.size(); ++i) {
    double v;
    m.evaluate( pd, handles[i], v, err );
    CPPUNIT_ASSERT(!err);
    ev1 += v;
    ev2 += v*v;
  }
  ev1 /= handles.size();
  ev2 /= handles.size();
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev1, v1, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev2, v2, 1e-6 );
}


void PMeanPMetricTest::test_indices() 
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m2( 2.0, &m );
  
  double v1, v2;
  msq_std::vector<size_t> indices;
  m2.evaluate_with_indices( pd, 0, v1, indices, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate( pd, 0, v2, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2, v1, 1e-6 );
  
  msq_std::vector<size_t> vertices;
  pd.element_by_index(0).get_vertex_indices( vertices );
  
  msq_std::sort( vertices.begin(), vertices.end() );
  msq_std::sort( indices.begin(), indices.end() );
  CPPUNIT_ASSERT( vertices == indices );
}

void PMeanPMetricTest::test_gradient()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m1( 1.0, &m );
  ElementPMeanP m2( 2.0, &m );
  
    // get vertices for later
  msq_std::vector<size_t> vertices;
  pd.element_by_index(0).get_vertex_indices( vertices );
  
    // evaluate without gradients
  double v1, v2, v3, v4;
  msq_std::vector<size_t> indices1, indices2, indices3, indices4;
  msq_std::vector<Vector3D> grads1, grads2;
  m1.evaluate_with_indices( pd, 0, v1, indices1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate_with_indices( pd, 0, v2, indices2, err );
  CPPUNIT_ASSERT(!err);
  
    // evaluate with gradients
  m1.evaluate_with_gradient( pd, 0, v3, indices3, grads1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate_with_gradient( pd, 0, v4, indices4, grads2, err );
  CPPUNIT_ASSERT(!err);
  
    // compare value and indices to eval w/out gradient
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v1, v3, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2, v4, 1e-6 );
  msq_std::sort( indices1.begin(), indices1.end() );
  msq_std::sort( indices3.begin(), indices3.end() );
  CPPUNIT_ASSERT( indices3 == indices1 );
  msq_std::sort( indices2.begin(), indices2.end() );
  msq_std::sort( indices4.begin(), indices4.end() );
  CPPUNIT_ASSERT( indices4 == indices2 );
  
    // setup evaluation of underying metric
  msq_std::vector<size_t> handles;
  m.get_element_evaluations( pd, 0, handles, err );
  CPPUNIT_ASSERT(!err);
  
    // calculate expected gradients
  msq_std::vector<size_t>::iterator j;
  msq_std::vector<Vector3D> expected1, expected2, temp;
  expected1.resize( vertices.size(), Vector3D(0,0,0) );
  expected2.resize( vertices.size(), Vector3D(0,0,0) );
  for (unsigned i = 0; i < handles.size(); ++i) {
    double v;
    m.evaluate_with_gradient( pd, handles[i], v, indices3, temp, err );
    CPPUNIT_ASSERT(!err);
    for (unsigned i = 0; i < indices3.size(); ++i) {
      j = std::find( vertices.begin(), vertices.end(), indices3[i] );
      CPPUNIT_ASSERT( j != vertices.end() );
      unsigned k = j - vertices.begin();
      expected1[k] += temp[i];
      expected2[k] += 2 * v * temp[i];
    }
  }
  for (unsigned i = 0; i < vertices.size(); ++i) {
    expected1[i] /= handles.size();
    expected2[i] /= handles.size();
  }
  
    // compare gradients
  for (unsigned i = 0; i < indices1.size(); ++i) {
    j = std::find( vertices.begin(), vertices.end(), indices1[i] );
    CPPUNIT_ASSERT( j != vertices.end() );
    unsigned k = j - vertices.begin();
    CPPUNIT_ASSERT_VECTORS_EQUAL( expected1[k], grads1[i], 1e-6 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( expected2[k], grads2[i], 1e-6 );
  }
}
