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
#else
# include <iostream>
#endif

using namespace Mesquite;
using namespace std;
const double epsilon = 1e-6;
#define ASSERT_VALUES_EQUAL( v1, v2, location, bits ) \
  ASSERT_MESSAGE( value_message( (location), (bits), (v1), (v2) ), \
                          (fabs((v1) - (v2)) < epsilon) )

inline const char* bintostr( unsigned bits )
{
  static char buffer_mem[6];
  char* buffer = buffer_mem;
  for (int i = sizeof(buffer_mem)-2; i >= 0; --i) {
    *buffer = (bits & (1<<i)) ? '1' : '0';
    ++buffer;
  }
  *buffer='\0';
  return buffer_mem;
}

inline CppUnit::Message value_message( unsigned location, unsigned bits, double v1, double v2 )
{
  char buffer[128];
  CppUnit::Message m( "equality assertion failed" );

  sprintf(buffer, "Expected : %f", v1 );
  m.addDetail( buffer );

  sprintf(buffer, "Actual   : %f", v2 );
  m.addDetail( buffer );

  if (location < 4) 
    sprintf(buffer, "Location : Corner %u", location );
  else if (location < 8)
    sprintf(buffer, "Location : Edge %u", location-4 );
  else if (location == 8)
    sprintf(buffer, "Location : Mid-element" );
  else
    sprintf(buffer, "Location : INVALID!!" );
  m.addDetail( buffer );

  sprintf(buffer, "Node Bits: %s", bintostr(bits) );
  m.addDetail( buffer );
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
    
    void test_corner_coeff( int corner, unsigned nodebits );
    void test_edge_coeff( int edge, unsigned nodebits );
    void test_mid_coeff( unsigned nodebits );
    
    void test_corner_derivs( int corner, unsigned nodebits );
    void test_edge_derivs( int edge, unsigned nodebits );
    void test_mid_derivs( unsigned nodebits );
    
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
                  unsigned nodebits, 
                  double xi, double eta, 
                  double results[9] )
{
    // Initial values are a) linear values for corners,
    // b) values for mid edges assuming no mid-face node,
    // and c) the mid-face value for a fully quadratic element.
  for (unsigned i = 0; i < 9; ++i)
    results[i] = function( i, xi, eta );
  
    // if center node is present, adjust mid-edge coefficients
  if (!(nodebits & 16u)) {
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
    if (!(nodebits & (1<<i))) {
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

static void get_coeffs( unsigned nodebits, double xi, double eta, 
                        double coeffs_out[9] )
  { eval( &N, nodebits, xi, eta, coeffs_out ); }

static void get_partial_wrt_xi( unsigned nodebits, double xi, double eta, 
                                double derivs_out[9] )
  { eval( &dNdxi, nodebits, xi, eta, derivs_out ); }

static void get_partial_wrt_eta( unsigned nodebits, double xi, double eta, 
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

static void check_valid_indices( std::vector<size_t> vertices, unsigned bits )
{
    // check valid size of list (at least three, at most all nodes)
  CPPUNIT_ASSERT( vertices.size() <= 9 );
  CPPUNIT_ASSERT( vertices.size() >= 3 );
    // make sure vertex indices are valid (in [0,8])
  std::sort( vertices.begin(), vertices.end() );
  CPPUNIT_ASSERT( vertices.back() <= 8 ); // max value less than 9
    // make sure there are no duplicates in the list
  std::vector<size_t>::iterator iter;
  iter = msq_std::unique( vertices.begin(), vertices.end() );
  CPPUNIT_ASSERT( iter == vertices.end() );

    // make all vertices are present in element
  for (unsigned i = 0; i < vertices.size(); ++i)
    if (vertices[i] > 4)
      CPPUNIT_ASSERT( bits & (1<<(vertices[i]-4)) );
}

static void check_no_zeros( std::vector<double> derivs )
{
  CPPUNIT_ASSERT_EQUAL( (size_t)0, derivs.size()%2 );
  for (unsigned i = 0; i < derivs.size()/2; ++i) {
    double dxi = derivs[2*i];
    double deta = derivs[2*i+1];
    CPPUNIT_ASSERT( fabs(dxi) > 1e-6 || fabs(deta) > 1e-6 );
  }
}

static void compare_coefficients( std::vector<double> coeffs,
                           const double* expected_coeffs,
                           unsigned loc, unsigned bits )
{
    // if vertex is not present, it must have a zero coefficient
  CPPUNIT_ASSERT_EQUAL( (size_t)9, coeffs.size() );
  if (!(bits & 1))
    ASSERT_VALUES_EQUAL( 0.0, coeffs[4], loc, bits );
  if (!(bits & 2))
    ASSERT_VALUES_EQUAL( 0.0, coeffs[5], loc, bits );
  if (!(bits & 4))
    ASSERT_VALUES_EQUAL( 0.0, coeffs[6], loc, bits );
  if (!(bits & 8))
    ASSERT_VALUES_EQUAL( 0.0, coeffs[7], loc, bits );
  if (!(bits & 16))
    ASSERT_VALUES_EQUAL( 0.0, coeffs[8], loc, bits );
    
    // compare expected and actual coefficient values
  CPPUNIT_ASSERT_EQUAL( (size_t)9, coeffs.size() );
  ASSERT_VALUES_EQUAL( expected_coeffs[0], coeffs[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[1], coeffs[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[2], coeffs[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[3], coeffs[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[4], coeffs[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[5], coeffs[5], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[6], coeffs[6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[7], coeffs[7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[8], coeffs[8], loc, bits );
}

static void compare_derivatives( std::vector<size_t> vertices,
                                 std::vector<double> derivs,
                                 const double* expected_dxi,
                                 const double* expected_deta,
                                 unsigned loc, unsigned bits )
{
  CPPUNIT_ASSERT_EQUAL( 2*vertices.size(), derivs.size() );
  check_valid_indices( vertices, bits );
  check_no_zeros( derivs );
  
    // Input has values in dxi & deta only for nodes in 'vertices'
    // Convert to values for every possible node, with zero's for
    // nodes that are not present.
  std::vector<double> expanded_dxi(9, 0.0), expanded_deta(9, 0.0);
  for (unsigned i = 0; i < vertices.size(); ++i) {
    expanded_dxi [vertices[i]] = derivs[2*i  ];
    expanded_deta[vertices[i]] = derivs[2*i+1];
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

void QuadLagrangeShapeTest::test_corner_coeff( int corner, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[9];
  get_coeffs( nodebits, corners[corner][XI], corners[corner][ETA], expected );
  
  std::vector<double> coeff;
  sf.coefficients_at_corner( corner, nodebits, coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, expected, corner, nodebits );
}

void QuadLagrangeShapeTest::test_edge_coeff( int edge, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[9];
  get_coeffs( nodebits, midedge[edge][XI], midedge[edge][ETA], expected );
  
  std::vector<double> coeff;
  sf.coefficients_at_mid_edge( edge, nodebits, coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, expected, edge+4, nodebits );
}

void QuadLagrangeShapeTest::test_mid_coeff( unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[9];
  get_coeffs( nodebits, midelem[XI], midelem[ETA], expected );
  
  std::vector<double> coeff;
  sf.coefficients_at_mid_elem( nodebits, coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, expected, 8, nodebits );
}

void QuadLagrangeShapeTest::test_corner_derivs( int corner, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected_dxi[9], expected_deta[9];
  get_partial_wrt_xi ( nodebits, corners[corner][XI], corners[corner][ETA], expected_dxi );
  get_partial_wrt_eta( nodebits, corners[corner][XI], corners[corner][ETA], expected_deta);
  
  std::vector<size_t> vertices;
  std::vector<double> derivs;
  sf.derivatives_at_corner( corner, nodebits, vertices, derivs, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, derivs, expected_dxi, expected_deta, corner, nodebits );
}

void QuadLagrangeShapeTest::test_edge_derivs( int edge, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected_dxi[9], expected_deta[9];
  get_partial_wrt_xi ( nodebits, midedge[edge][XI], midedge[edge][ETA], expected_dxi );
  get_partial_wrt_eta( nodebits, midedge[edge][XI], midedge[edge][ETA], expected_deta);
  
  std::vector<size_t> vertices;
  std::vector<double> derivs;
  sf.derivatives_at_mid_edge( edge, nodebits, vertices, derivs, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, derivs, expected_dxi, expected_deta, edge+4, nodebits );
}

void QuadLagrangeShapeTest::test_mid_derivs( unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected_dxi[9], expected_deta[9];
  get_partial_wrt_xi ( nodebits, midelem[XI], midelem[ETA], expected_dxi );
  get_partial_wrt_eta( nodebits, midelem[XI], midelem[ETA], expected_deta);
  
  std::vector<size_t> vertices;
  std::vector<double> derivs;
  sf.derivatives_at_mid_elem( nodebits, vertices, derivs, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, derivs, expected_dxi, expected_deta, 8, nodebits );
}

void QuadLagrangeShapeTest::test_coeff_corners()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every corner
    for (unsigned i = 0; i < 4; ++i) 
      test_corner_coeff( i, j );
}

void QuadLagrangeShapeTest::test_coeff_edges()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every edge
    for (unsigned i = 0; i < 4; ++i) 
      test_edge_coeff( i, j );
}

void QuadLagrangeShapeTest::test_coeff_center()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
    test_mid_coeff( j );
}

void QuadLagrangeShapeTest::test_deriv_corners()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every corner
    for (unsigned i = 0; i < 4; ++i) 
      test_corner_derivs( i, j );
}

void QuadLagrangeShapeTest::test_deriv_edges()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
      // for every edge
    for (unsigned i = 0; i < 4; ++i) 
      test_edge_derivs( i, j );
}

void QuadLagrangeShapeTest::test_deriv_center()
{
    // for every possible combination of higher-order nodes
    // (0x1F = 11111 : five possible higher-order nodes in quad)
  for (unsigned j = 0; j <= 0x1Fu; ++j)
    test_mid_derivs( j );
}
