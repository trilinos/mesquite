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


/** \file QuadLagrangeShapeTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "QuadLagrangeShape.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"

#include "UnitUtil.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <vector.h>
# include <algorithm.h>
#else
# include <vector>
# include <algorithm>
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
# include <iostream.h>
# include <sstream.h>
#else
# include <iostream>
# include <sstream>
#endif

using namespace Mesquite;
using namespace std;
const double epsilon = 1e-6;
#define ASSERT_VALUES_EQUAL( v1, v2, location, bits ) \
  ASSERT_MESSAGE( value_message( (location), (bits), (v1), (v2) ), \
                          (fabs((v1) - (v2)) < epsilon) )

static inline CppUnit::Message value_message( unsigned location, NodeSet bits, double v1, double v2 )
{
  CppUnit::Message m( "equality assertion failed" );

  msq_stdio::ostringstream buffer1;
  buffer1 << "Expected : " << v1;
   m.addDetail( buffer1.str() );

  msq_stdio::ostringstream buffer2;
  buffer2 << "Actual   : " << v2;
  m.addDetail( buffer2.str() );

  msq_stdio::ostringstream buffer3;
  buffer3 << "Location : ";
  if (location < 4) 
    buffer3 << "Corner " << location;
  else if (location < 8)
    buffer3 << "Edge " << location-4;
  else if (location == 8)
    buffer3 << "Mid-element";
  else
    buffer3 << "INVALID!!";
  m.addDetail( buffer3.str() );

  msq_stdio::ostringstream buffer4;
  buffer4 << "Node Bits: " << bits;
  m.addDetail( buffer4.str() );
  return m;
}

class QuadLagrangeShapeTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(QuadLagrangeShapeTest);

    CPPUNIT_TEST(test_coeff_corners);
    CPPUNIT_TEST(test_coeff_edges);
    CPPUNIT_TEST(test_coeff_center);

    CPPUNIT_TEST(test_deriv_corners);
    CPPUNIT_TEST(test_deriv_edges);
    CPPUNIT_TEST(test_deriv_center);
    
    CPPUNIT_TEST_SUITE_END();
  
    QuadLagrangeShape sf;
    
    void test_corner_coeff( int corner, NodeSet nodeset );
    void test_edge_coeff( int edge, NodeSet nodeset );
    void test_mid_coeff( NodeSet nodeset );
    
    void test_corner_derivs( int corner, NodeSet nodeset );
    void test_edge_derivs( int edge, NodeSet nodeset );
    void test_mid_derivs( NodeSet nodeset );
    
  public:

    void test_coeff_corners();
    void test_coeff_edges();
    void test_coeff_center();

    void test_deriv_corners();
    void test_deriv_edges();
    void test_deriv_center();
};



CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QuadLagrangeShapeTest, "QuadLagrangeShapeTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QuadLagrangeShapeTest, "Unit");

// 2nd order lagrange polynomial indices indexed by (vertex,xi/eta)
enum { XI = 0, ETA = 1 };
static double N0(double xi, double eta) { return 0.25*(1 - xi)*(1 - eta); }
static double N1(double xi, double eta) { return 0.25*(1 + xi)*(1 - eta); }
static double N2(double xi, double eta) { return 0.25*(1 + xi)*(1 + eta); }
static double N3(double xi, double eta) { return 0.25*(1 - xi)*(1 + eta); }
static double N4(double xi, double eta) { return 0.5*(1 - xi*xi)*(1 - eta); }
static double N5(double xi, double eta) { return 0.5*(1 - eta*eta)*(1 + xi); }
static double N6(double xi, double eta) { return 0.5*(1 - xi*xi)*(1 + eta); }
static double N7(double xi, double eta) { return 0.5*(1 - eta*eta)*(1 - xi); }
static double N8(double xi, double eta) { return (1 - xi*xi)*(1 - eta*eta); }

static double dN0dxi(double xi, double eta) { return -0.25*(1 - eta); }
static double dN1dxi(double xi, double eta) { return  0.25*(1 - eta); }
static double dN2dxi(double xi, double eta) { return  0.25*(1 + eta); }
static double dN3dxi(double xi, double eta) { return -0.25*(1 + eta); }
static double dN4dxi(double xi, double eta) { return -xi*(1 - eta); }
static double dN5dxi(double xi, double eta) { return  0.5*(1 - eta*eta); }
static double dN6dxi(double xi, double eta) { return -xi*(1 + eta); }
static double dN7dxi(double xi, double eta) { return -0.5*(1 - eta*eta); }
static double dN8dxi(double xi, double eta) { return -2.0*xi*(1 - eta*eta); }

static double dN0deta(double xi, double eta) { return -0.25*(1 - xi); }
static double dN1deta(double xi, double eta) { return -0.25*(1 + xi); }
static double dN2deta(double xi, double eta) { return  0.25*(1 + xi); }
static double dN3deta(double xi, double eta) { return  0.25*(1 - xi); }
static double dN4deta(double xi, double eta) { return -0.5*(1 - xi*xi); }
static double dN5deta(double xi, double eta) { return -eta*(1 + xi); }
static double dN6deta(double xi, double eta) { return  0.5*(1 - xi*xi); }
static double dN7deta(double xi, double eta) { return -eta*(1 - xi); }
static double dN8deta(double xi, double eta) { return -2.0*eta*(1 - xi*xi); }

typedef double (*N_t)(double, double);
static const N_t N_a[9] = { &N0, &N1, &N2, &N3, &N4, &N5, &N6, &N7, &N8 };
static const N_t dNdxi_a[9] = { &dN0dxi, &dN1dxi, &dN2dxi, &dN3dxi, &dN4dxi, &dN5dxi, &dN6dxi, &dN7dxi, &dN8dxi };
static const N_t dNdeta_a[9] = { &dN0deta, &dN1deta, &dN2deta, &dN3deta, &dN4deta, &dN5deta, &dN6deta, &dN7deta, &dN8deta };

// Mapping function coefficient functions (i<-vertex number)
static double N( unsigned i, double xi, double eta )
  { return N_a[i](xi,eta); }

// Partial derivatives of mapping function coefficients
static double dNdxi( unsigned i, double xi, double eta )
  { return dNdxi_a[i](xi,eta); }
static double dNdeta( unsigned i, double xi, double eta )
  { return dNdeta_a[i](xi,eta); }


// Evaluate one of N, dNdxi or dNdeta for each vertex
typedef double (*f_t)(unsigned, double, double);
static void eval( f_t function, 
                  NodeSet nodeset, 
                  double xi, double eta, 
                  double results[9] )
{
    // Initial values are a) linear values for corners,
    // b) values for mid edges assuming no mid-face node,
    // and c) the mid-face value for a fully quadratic element.
  for (unsigned i = 0; i < 9; ++i)
    results[i] = function( i, xi, eta );
  
    // if center node is present, adjust mid-edge coefficients
  if (!nodeset.mid_face_node(0)) {
    results[8] = 0;
  }
  else {
    for (unsigned i = 0; i < 4; ++i)
      results[i] -= 0.25 * results[8];
    for (unsigned i = 4; i < 8; ++i)
      results[i] -= 0.5 * results[8];
  }
  
    // if mid-edge nodes are present, adjust values for adjacent corners
  for (unsigned i = 0; i < 4; ++i) {
    if (!nodeset.mid_edge_node(i)) {
      results[i+4] = 0.0;
    }
    else {
      results[ i     ] -= 0.5 * results[i+4]; // 1st adjacent corner
      results[(i+1)%4] -= 0.5 * results[i+4]; // 2nd adjacent corner
    }
  }
}

// Finally, what all the above stuff was building up to:
// functions to query mapping function.

static void get_coeffs( NodeSet nodebits, double xi, double eta, 
                        double coeffs_out[9] )
  { eval( &N, nodebits, xi, eta, coeffs_out ); }

static void get_partial_wrt_xi( NodeSet nodebits, double xi, double eta, 
                                double derivs_out[9] )
  { eval( &dNdxi, nodebits, xi, eta, derivs_out ); }

static void get_partial_wrt_eta( NodeSet nodebits, double xi, double eta, 
                                double derivs_out[9] )
  { eval( &dNdeta, nodebits, xi, eta, derivs_out ); }


// Pre-defined sample points (xi and eta values)
static const double corners[4][2] = { { -1, -1 },
                                      {  1, -1 },
                                      {  1,  1 },
                                      { -1,  1 } };
static const double midedge[4][2] = { {  0, -1 },
                                      {  1,  0 },
                                      {  0,  1 },
                                      { -1,  0 } };
static const double midelem[2] = { 0, 0 };

static void check_valid_indices( const size_t* vertices, size_t num_vtx, NodeSet nodeset )
{
    // check valid size of list (at least three, at most all nodes)
  CPPUNIT_ASSERT( num_vtx <= 9 );
  CPPUNIT_ASSERT( num_vtx >= 3 );
    // make sure vertex indices are valid (in [0,8])
  size_t vertcopy[9];
  std::copy( vertices, vertices+num_vtx, vertcopy );
  std::sort( vertcopy, vertcopy+num_vtx );
  CPPUNIT_ASSERT( vertcopy[num_vtx-1] <= 8 ); // max value less than 9
    // make sure there are no duplicates in the list
  const size_t* iter = msq_std::unique( vertcopy, vertcopy+num_vtx );
  CPPUNIT_ASSERT( iter == vertcopy+num_vtx );

    // make all vertices are present in element
  for (unsigned i = 0; i < num_vtx; ++i) {
    if (vertcopy[i] == 8)
      CPPUNIT_ASSERT( nodeset.mid_face_node(0) );
    else if (vertcopy[i] >= 4)
      CPPUNIT_ASSERT( nodeset.mid_edge_node( vertcopy[i] - 4 ) );
  }
}

static void check_no_zeros( const MsqVector<2>* derivs, size_t num_vtx )
{
  for (unsigned i = 0; i < num_vtx; ++i) {
    double dxi = derivs[i][0];
    double deta = derivs[i][1];
    CPPUNIT_ASSERT( fabs(dxi) > 1e-6 || fabs(deta) > 1e-6 );
  }
}

static void compare_coefficients( const double* coeffs,
                                  const size_t* indices,
                                  size_t num_coeff,
                                  const double* expected_coeffs,
                                  unsigned loc, NodeSet bits )
{
    // find the location in the returned list for each node
  size_t revidx[9];
  double test_vals[9];
  for (size_t i = 0; i < 9; ++i) {
    revidx[i] = std::find( indices, indices+num_coeff, i ) - indices;
    test_vals[i] = (revidx[i] == num_coeff) ? 0.0 : coeffs[revidx[i]];
  }

    // Check that index list doesn't contain any nodes not actually
    // present in the element.
  CPPUNIT_ASSERT( bits.mid_edge_node(0) || (revidx[4] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(1) || (revidx[5] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(2) || (revidx[6] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(3) || (revidx[7] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_face_node(0) || (revidx[8] == num_coeff) );
    
    // compare expected and actual coefficient values
  ASSERT_VALUES_EQUAL( expected_coeffs[0], test_vals[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[1], test_vals[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[2], test_vals[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[3], test_vals[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[4], test_vals[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[5], test_vals[5], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[6], test_vals[6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[7], test_vals[7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[8], test_vals[8], loc, bits );
}

static void compare_derivatives( const size_t* vertices,
                                 size_t num_vtx,
                                 const MsqVector<2>* derivs,
                                 const double* expected_dxi,
                                 const double* expected_deta,
                                 unsigned loc, NodeSet bits )
{
  check_valid_indices( vertices, num_vtx, bits );
  check_no_zeros( derivs, num_vtx );
  
    // Input has values in dxi & deta only for nodes in 'vertices'
    // Convert to values for every possible node, with zero's for
    // nodes that are not present.
  double expanded_dxi[9], expanded_deta[9];
  std::fill( expanded_dxi, expanded_dxi+9, 0.0 );
  std::fill( expanded_deta, expanded_deta+9, 0.0 );
  for (unsigned i = 0; i < num_vtx; ++i) {
    expanded_dxi [vertices[i]] = derivs[i][0];
    expanded_deta[vertices[i]] = derivs[i][1];
  }
  
  ASSERT_VALUES_EQUAL( expected_dxi[0], expanded_dxi[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[1], expanded_dxi[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[2], expanded_dxi[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[3], expanded_dxi[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[4], expanded_dxi[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[5], expanded_dxi[5], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[6], expanded_dxi[6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[7], expanded_dxi[7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_dxi[8], expanded_dxi[8], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_deta[0], expanded_deta[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[1], expanded_deta[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[2], expanded_deta[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[3], expanded_deta[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[4], expanded_deta[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[5], expanded_deta[5], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[6], expanded_deta[6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[7], expanded_deta[7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_deta[8], expanded_deta[8], loc, bits );
}

void QuadLagrangeShapeTest::test_corner_coeff( int corner, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[9];
  get_coeffs( nodebits, corners[corner][XI], corners[corner][ETA], expected );
  
  double coeff[100];
  size_t num_coeff = 11, indices[100];
  sf.coefficients( Sample(0, corner), nodebits, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, num_coeff, expected, corner, nodebits );
}

void QuadLagrangeShapeTest::test_edge_coeff( int edge, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[9];
  get_coeffs( nodebits, midedge[edge][XI], midedge[edge][ETA], expected );
  
  double coeff[100];
  size_t num_coeff = 11, indices[100];
  sf.coefficients( Sample(1, edge), nodebits, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, num_coeff, expected, edge+4, nodebits );
}

void QuadLagrangeShapeTest::test_mid_coeff( NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[9];
  get_coeffs( nodebits, midelem[XI], midelem[ETA], expected );
  
  double coeff[100];
  size_t num_coeff = 11, indices[100];
  sf.coefficients( Sample(2, 0), nodebits, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, num_coeff, expected, 8, nodebits );
}

void QuadLagrangeShapeTest::test_corner_derivs( int corner, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected_dxi[9], expected_deta[9];
  get_partial_wrt_xi ( nodebits, corners[corner][XI], corners[corner][ETA], expected_dxi );
  get_partial_wrt_eta( nodebits, corners[corner][XI], corners[corner][ETA], expected_deta);
  
  size_t vertices[100], num_vtx = 23;
  MsqVector<2> derivs[100];
  sf.derivatives( Sample(0, corner), nodebits, vertices, derivs, num_vtx, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, num_vtx, derivs, expected_dxi, expected_deta, corner, nodebits );
}

void QuadLagrangeShapeTest::test_edge_derivs( int edge, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected_dxi[9], expected_deta[9];
  get_partial_wrt_xi ( nodebits, midedge[edge][XI], midedge[edge][ETA], expected_dxi );
  get_partial_wrt_eta( nodebits, midedge[edge][XI], midedge[edge][ETA], expected_deta);
  
  size_t vertices[100], num_vtx = 23;
  MsqVector<2> derivs[100];
  sf.derivatives( Sample(1, edge), nodebits, vertices, derivs, num_vtx, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, num_vtx, derivs, expected_dxi, expected_deta, edge+4, nodebits );
}

void QuadLagrangeShapeTest::test_mid_derivs( NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected_dxi[9], expected_deta[9];
  get_partial_wrt_xi ( nodebits, midelem[XI], midelem[ETA], expected_dxi );
  get_partial_wrt_eta( nodebits, midelem[XI], midelem[ETA], expected_deta);
  
  size_t vertices[100], num_vtx = 23;
  MsqVector<2> derivs[100];
  sf.derivatives( Sample(2, 0), nodebits, vertices, derivs, num_vtx, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, num_vtx, derivs, expected_dxi, expected_deta, 8, nodebits );
}

static NodeSet nodeset_from_bits( unsigned bits )
{
  NodeSet result;
  for (unsigned i = 0; i < 4; ++i)
    if (bits & (1<<i))
      result.set_mid_edge_node(i);
  if (bits & (1<<4))
    result.set_mid_face_node(0);
  return result;
}

void QuadLagrangeShapeTest::test_coeff_corners()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every corner
    for (unsigned i = 0; i < 4; ++i) 
      test_corner_coeff( i, nodeset_from_bits(j) );
}

void QuadLagrangeShapeTest::test_coeff_edges()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every edge
    for (unsigned i = 0; i < 4; ++i) 
      test_edge_coeff( i, nodeset_from_bits(j) );
}

void QuadLagrangeShapeTest::test_coeff_center()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
    test_mid_coeff( nodeset_from_bits(j) );
}

void QuadLagrangeShapeTest::test_deriv_corners()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every corner
    for (unsigned i = 0; i < 4; ++i) 
      test_corner_derivs( i, nodeset_from_bits(j));
}

void QuadLagrangeShapeTest::test_deriv_edges()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every edge
    for (unsigned i = 0; i < 4; ++i) 
      test_edge_derivs( i, nodeset_from_bits(j) );
}

void QuadLagrangeShapeTest::test_deriv_center()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
    test_mid_derivs( nodeset_from_bits(j) );
}
