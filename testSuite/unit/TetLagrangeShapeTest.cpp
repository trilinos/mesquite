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


/** \file TetLagrangeShapeTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TetLagrangeShape.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"
#include <math.h>

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
  static char buffer_mem[12];
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
  else if (location < 10)
    sprintf(buffer, "Location : Edge %u", location-4 );
  else if (location < 14)
    sprintf(buffer, "Location : Face %u", location-10 );
  else if (location == 14)
    sprintf(buffer, "Location : Mid-element" );
  else
    sprintf(buffer, "Invalid location" );
  m.addDetail( buffer );

  sprintf(buffer, "Node Bits: %s", bintostr(bits) );
  m.addDetail( buffer );
  return m;
}

class TetLagrangeShapeTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(TetLagrangeShapeTest);

    CPPUNIT_TEST(test_coeff_corners);
    CPPUNIT_TEST(test_coeff_edges);
    CPPUNIT_TEST(test_coeff_faces);
    CPPUNIT_TEST(test_coeff_center);

    CPPUNIT_TEST(test_deriv_corners);
    CPPUNIT_TEST(test_deriv_edges);
    CPPUNIT_TEST(test_deriv_faces);
    CPPUNIT_TEST(test_deriv_center);

    CPPUNIT_TEST(test_mid_elem_node_coeff);
    CPPUNIT_TEST(test_mid_elem_node_deriv);
    CPPUNIT_TEST(test_mid_face_node_coeff);
    CPPUNIT_TEST(test_mid_face_node_deriv);
    
    CPPUNIT_TEST_SUITE_END();
  
    TetLagrangeShape sf;
    
    void test_corner_coeff( int corner, unsigned nodebits );
    void test_edge_coeff( int edge, unsigned nodebits );
    void test_face_coeff( int face, unsigned nodebits );
    void test_mid_coeff( unsigned nodebits );
    
    void test_corner_derivs( int corner, unsigned nodebits );
    void test_edge_derivs( int edge, unsigned nodebits );
    void test_face_derivs( int face, unsigned nodebits );
    void test_mid_derivs( unsigned nodebits );
    
    void test_invalid_nodebits_coeff( unsigned nodebits );
    void test_invalid_nodebits_deriv( unsigned nodebits );
    
  public:

    void test_coeff_corners();
    void test_coeff_edges();
    void test_coeff_faces();
    void test_coeff_center();

    void test_deriv_corners();
    void test_deriv_edges();
    void test_deriv_faces();
    void test_deriv_center();
    
    void test_mid_elem_node_coeff();
    void test_mid_elem_node_deriv();
    void test_mid_face_node_coeff();
    void test_mid_face_node_deriv();
};



CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TetLagrangeShapeTest, "TetLagrangeShapeTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TetLagrangeShapeTest, "Unit");

static double N0( double r, double  , double   ) { return r*(2*r - 1); }
static double N1( double  , double s, double   ) { return s*(2*s - 1); }
static double N2( double  , double  , double t ) { return t*(2*t - 1); }
static double N3( double r, double s, double t ) { double u = 1 - r - s - t; return u*(2*u - 1); }
static double N4( double r, double s, double   ) { return 4*r*s; }
static double N5( double  , double s, double t ) { return 4*s*t; }
static double N6( double r, double  , double t ) { return 4*r*t; }
static double N7( double r, double s, double t ) { double u = 1 - r - s - t; return 4*r*u; }
static double N8( double r, double s, double t ) { double u = 1 - r - s - t; return 4*s*u; }
static double N9( double r, double s, double t ) { double u = 1 - r - s - t; return 4*t*u; }

static double dN0dr( double r, double  , double   ) { return 4*r - 1; }
static double dN0ds( double  , double  , double   ) { return 0; }
static double dN0dt( double  , double  , double   ) { return 0; }

static double dN1dr( double  , double  , double   ) { return 0; }
static double dN1ds( double  , double s, double   ) { return 4*s - 1; }
static double dN1dt( double  , double  , double   ) { return 0; }

static double dN2dr( double  , double  , double   ) { return 0; }
static double dN2ds( double  , double  , double   ) { return 0; }
static double dN2dt( double  , double  , double t ) { return 4*t - 1; }

static double dN3dr( double r, double s, double t ) { return 4*(r + s + t) - 3; }
static double dN3ds( double r, double s, double t ) { return 4*(r + s + t) - 3; }
static double dN3dt( double r, double s, double t ) { return 4*(r + s + t) - 3; }

static double dN4dr( double  , double s, double   ) { return 4*s; }
static double dN4ds( double r, double  , double   ) { return 4*r; }
static double dN4dt( double  , double  , double   ) { return 0; }

static double dN5dr( double  , double  , double   ) { return 0; }
static double dN5ds( double  , double  , double t ) { return 4*t; }
static double dN5dt( double  , double s, double   ) { return 4*s; }

static double dN6dr( double  , double  , double t ) { return 4*t; }
static double dN6ds( double  , double  , double   ) { return 0; }
static double dN6dt( double r, double  , double   ) { return 4*r; }

static double dN7dr( double r, double s, double t ) { return 4*(1 - 2*r - s - t); }
static double dN7ds( double r, double  , double   ) { return -4*r; }
static double dN7dt( double r, double  , double   ) { return -4*r; }

static double dN8dr( double  , double s, double   ) { return -4*s; }
static double dN8ds( double r, double s, double t ) { return 4*(1 - r - 2*s - t); }
static double dN8dt( double  , double s, double   ) { return -4*s; }

static double dN9dr( double  , double  , double t ) { return -4*t; }
static double dN9ds( double  , double  , double t ) { return -4*t; }
static double dN9dt( double r, double s, double t ) { return 4*(1 - r - s - 2*t); }

typedef double (*N_t)(double,double,double);
static const N_t N[] = { &N0, &N1, &N2, &N3, &N4, &N5, &N6, &N7, &N8, &N9 };
static const N_t dNdr[] = { &dN0dr, &dN1dr, &dN2dr, &dN3dr, &dN4dr, 
                     &dN5dr, &dN6dr, &dN7dr, &dN8dr, &dN9dr };
static const N_t dNds[] = { &dN0ds, &dN1ds, &dN2ds, &dN3ds, &dN4ds, 
                     &dN5ds, &dN6ds, &dN7ds, &dN8ds, &dN9ds };
static const N_t dNdt[] = { &dN0dt, &dN1dt, &dN2dt, &dN3dt, &dN4dt, 
                     &dN5dt, &dN6dt, &dN7dt, &dN8dt, &dN9dt };

static const double rst_corner[][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
static const double rst_edge[][3] = { {0.5, 0.5, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5},
                               {0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5} };
static const double rst_face[][3] = { {1./3, 1./3, 0.00}, {0.00, 1./3, 1./3}, 
                               {1./3, 0.00, 1./3}, {1./3, 1./3, 1./3} };
static const double rst_mid[3] = { 0.25, 0.25, 0.25 };

static unsigned edges[][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 },
                        { 0, 3 }, { 1, 3 }, { 2, 3 } };

static void get_coeff( unsigned nodebits, const double* rst, double* coeffs )
{
  for (int i = 0; i < 10; ++i) 
    coeffs[i] = (*N[i])(rst[0], rst[1], rst[2]);
  for (int i = 0; i < 6; ++i) 
    if (!(nodebits & 1<<i)) {
      coeffs[edges[i][0]] += 0.5 * coeffs[i+4];
      coeffs[edges[i][1]] += 0.5 * coeffs[i+4];
      coeffs[i+4] = 0;
    }
}

static void get_derivs( unsigned nodebits, const double* rst, double* derivs )
{
  for (int i = 0; i < 10; ++i) {
    derivs[3*i  ] = (*dNdr[i])(rst[0], rst[1], rst[2]);
    derivs[3*i+1] = (*dNds[i])(rst[0], rst[1], rst[2]);
    derivs[3*i+2] = (*dNdt[i])(rst[0], rst[1], rst[2]);
  }
  for (int i = 0; i < 6; ++i) 
    if (!(nodebits & 1<<i)) {
      int j = edges[i][0];
      derivs[3*j  ] += 0.5 * derivs[3*i+12];
      derivs[3*j+1] += 0.5 * derivs[3*i+13];
      derivs[3*j+2] += 0.5 * derivs[3*i+14];
      j = edges[i][1];
      derivs[3*j  ] += 0.5 * derivs[3*i+12];
      derivs[3*j+1] += 0.5 * derivs[3*i+13];
      derivs[3*j+2] += 0.5 * derivs[3*i+14];
      derivs[3*i+12] = 0.0;
      derivs[3*i+13] = 0.0;
      derivs[3*i+14] = 0.0;
    }
}

static void check_valid_indices( const size_t* vertices, size_t num_vtx )
{
  CPPUNIT_ASSERT( num_vtx < 11 );
  CPPUNIT_ASSERT( num_vtx > 3 );
  size_t vtxcopy[10];
  std::copy( vertices, vertices+num_vtx, vtxcopy );
  std::sort( vtxcopy, vtxcopy+num_vtx );
  for (unsigned i = 1; i < num_vtx; ++i) {
    CPPUNIT_ASSERT( vtxcopy[i] != vtxcopy[i-1] );
    CPPUNIT_ASSERT( vtxcopy[i] < 10 );
  }
}

static void check_no_zeros( const double* derivs, size_t num_vtx )
{
  unsigned i = 0;
  while (i < 3*num_vtx) {
    double dr = derivs[i++]; 
    double ds = derivs[i++]; 
    double dt = derivs[i++]; 
    CPPUNIT_ASSERT( (fabs(dr) > 1e-6) || (fabs(ds) > 1e-6) || (fabs(dt) > 1e-6) );
  }
}

static void compare_coefficients( const double* coeffs,
                                  const double* expected_coeffs,
                                  size_t num_coeff,
                                  unsigned loc, unsigned bits )
{
  CPPUNIT_ASSERT_EQUAL( (size_t)10, num_coeff );
  ASSERT_VALUES_EQUAL( expected_coeffs[0], coeffs[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[1], coeffs[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[2], coeffs[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[3], coeffs[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[4], coeffs[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[5], coeffs[5], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[6], coeffs[6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[7], coeffs[7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[8], coeffs[8], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[9], coeffs[9], loc, bits );
}

static void compare_derivatives( const size_t* vertices,
                                 size_t num_vtx,
                                 const double* derivs,
                                 const double* expected_derivs,
                                 unsigned loc, unsigned bits )
{
  check_valid_indices( vertices, num_vtx );
  check_no_zeros( derivs, num_vtx );
  double expanded_derivs[30];
  memset( expanded_derivs, 0, sizeof(expanded_derivs) );
  for (unsigned i = 0; i < num_vtx; ++i) {
    expanded_derivs[3*vertices[i]  ] = derivs[3*i  ];
    expanded_derivs[3*vertices[i]+1] = derivs[3*i+1];
    expanded_derivs[3*vertices[i]+2] = derivs[3*i+2];
  }
  
  ASSERT_VALUES_EQUAL( expected_derivs[ 0], expanded_derivs[ 0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 1], expanded_derivs[ 1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 2], expanded_derivs[ 2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[ 3], expanded_derivs[ 3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 4], expanded_derivs[ 4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 5], expanded_derivs[ 5], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[ 6], expanded_derivs[ 6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 7], expanded_derivs[ 7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 8], expanded_derivs[ 8], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[ 9], expanded_derivs[ 9], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[10], expanded_derivs[10], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[11], expanded_derivs[11], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[12], expanded_derivs[12], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[13], expanded_derivs[13], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[14], expanded_derivs[14], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[15], expanded_derivs[15], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[16], expanded_derivs[16], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[17], expanded_derivs[17], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[18], expanded_derivs[18], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[19], expanded_derivs[19], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[20], expanded_derivs[20], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[21], expanded_derivs[21], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[22], expanded_derivs[22], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[23], expanded_derivs[23], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[24], expanded_derivs[24], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[25], expanded_derivs[25], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[26], expanded_derivs[26], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[27], expanded_derivs[27], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[28], expanded_derivs[28], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[29], expanded_derivs[29], loc, bits );
}

void TetLagrangeShapeTest::test_corner_coeff( int corner, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_corner[corner], expected );
  
  double coeff[100];
  size_t n = 29;
  sf.coefficients_at_corner( corner, nodebits, coeff, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, expected, n, corner, nodebits );
}

void TetLagrangeShapeTest::test_edge_coeff( int edge, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_edge[edge], expected );
  
  double coeff[100];
  size_t n = 29;
  sf.coefficients_at_mid_edge( edge, nodebits, coeff, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, expected, n, edge+4, nodebits );
}

void TetLagrangeShapeTest::test_face_coeff( int face, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_face[face], expected );
  
  double coeff[100];
  size_t n = 29;
  sf.coefficients_at_mid_face( face, nodebits, coeff, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, expected, n, face+10, nodebits );
}

void TetLagrangeShapeTest::test_mid_coeff( unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_mid, expected );
  
  double coeff[100];
  size_t n = 29;
  sf.coefficients_at_mid_elem( nodebits, coeff, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, expected, n, 14, nodebits );
}

void TetLagrangeShapeTest::test_corner_derivs( int corner, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[30];
  get_derivs( nodebits, rst_corner[corner], expected );
  
  double derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives_at_corner( corner, nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, corner, nodebits );
}

void TetLagrangeShapeTest::test_edge_derivs( int edge, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[30];
  get_derivs( nodebits, rst_edge[edge], expected );
  
  double derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives_at_mid_edge( edge, nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, edge+4, nodebits );
}

void TetLagrangeShapeTest::test_face_derivs( int face, unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[30];
  get_derivs( nodebits, rst_face[face], expected );
  
  double derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives_at_mid_face( face, nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, face+10, nodebits );
}

void TetLagrangeShapeTest::test_mid_derivs( unsigned nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[30];
  get_derivs( nodebits, rst_mid, expected );
  
  double derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives_at_mid_elem( nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, 14, nodebits );
}

void TetLagrangeShapeTest::test_coeff_corners()
{
  for (unsigned i = 0; i < 0x3F; ++i) {
    test_corner_coeff( 0, i );
    test_corner_coeff( 1, i );
    test_corner_coeff( 2, i );
    test_corner_coeff( 3, i );
  }
}

void TetLagrangeShapeTest::test_coeff_edges()
{
  for (unsigned i = 0; i < 0x3F; ++i) {
    test_edge_coeff( 0, 0 );
    test_edge_coeff( 1, 0 );
    test_edge_coeff( 2, 0 );
    test_edge_coeff( 3, 0 );
    test_edge_coeff( 4, 0 );
    test_edge_coeff( 5, 0 );
  }
}

void TetLagrangeShapeTest::test_coeff_faces()
{
  for (unsigned i = 0; i < 0x3F; ++i) {
    test_face_coeff( 0, 0 );
    test_face_coeff( 1, 0 );
    test_face_coeff( 2, 0 );
    test_face_coeff( 3, 0 );
  }
}

void TetLagrangeShapeTest::test_coeff_center()
{
  for (unsigned i = 0; i < 0x3F; ++i) {
    test_mid_coeff( i );
  }
}

void TetLagrangeShapeTest::test_deriv_corners()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_corner_derivs( 0, i );
    test_corner_derivs( 1, i );
    test_corner_derivs( 2, i );
    test_corner_derivs( 3, i );
  }
}

void TetLagrangeShapeTest::test_deriv_edges()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_edge_derivs( 0, i );
    test_edge_derivs( 1, i );
    test_edge_derivs( 2, i );
    test_edge_derivs( 3, i );
    test_edge_derivs( 4, i );
    test_edge_derivs( 5, i );
  }
}

void TetLagrangeShapeTest::test_deriv_faces()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_face_derivs( 0, i );
    test_face_derivs( 1, i );
    test_face_derivs( 2, i );
    test_face_derivs( 3, i );
  }
}

void TetLagrangeShapeTest::test_deriv_center()
{
  for (unsigned i = 0; i < 0x3F; ++i) {
    test_mid_derivs( i );
  }
}

void TetLagrangeShapeTest::test_invalid_nodebits_coeff( unsigned bits )
{
  MsqError err;
  double coeff[100];
  size_t n;
  
  sf.coefficients_at_corner( 0, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_corner( 1, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_corner( 2, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_corner( 3, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.coefficients_at_mid_edge( 0, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_edge( 1, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_edge( 2, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_edge( 3, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_edge( 4, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_edge( 5, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.coefficients_at_mid_face( 0, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_face( 1, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_face( 2, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients_at_mid_face( 3, bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.coefficients_at_mid_elem( bits, coeff, n, err );
  CPPUNIT_ASSERT( err );
}

void TetLagrangeShapeTest::test_invalid_nodebits_deriv( unsigned bits )
{
  MsqError err;
  size_t verts[100], n;
  double derivs[100];
  
  sf.derivatives_at_corner( 0, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_corner( 1, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_corner( 2, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_corner( 3, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.derivatives_at_mid_edge( 0, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_edge( 1, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_edge( 2, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_edge( 3, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_edge( 4, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_edge( 5, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.derivatives_at_mid_face( 0, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_face( 1, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_face( 2, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives_at_mid_face( 3, bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.derivatives_at_mid_elem( bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
}


void TetLagrangeShapeTest::test_mid_elem_node_coeff()
{
  test_invalid_nodebits_coeff( 0x400 );
}

void TetLagrangeShapeTest::test_mid_elem_node_deriv()
{
  test_invalid_nodebits_deriv( 0x400 );
}

void TetLagrangeShapeTest::test_mid_face_node_coeff()
{
  test_invalid_nodebits_coeff( 0x040 );
  test_invalid_nodebits_coeff( 0x200 );
}

void TetLagrangeShapeTest::test_mid_face_node_deriv()
{
  test_invalid_nodebits_deriv( 0x040 );
  test_invalid_nodebits_deriv( 0x200 );
}

