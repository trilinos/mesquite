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


/** \file TriLagrangeShapeTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TriLagrangeShape.hpp"
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
  if (location < 3) 
    buffer3 << "Corner " << location;
  else if (location < 6)
    buffer3 << "Edge " << location-3;
  else if (location == 6)
    buffer3 << "Mid-element";
  else
    buffer3 << "INVALID!!";
  m.addDetail( buffer3.str() );

  msq_stdio::ostringstream buffer4;
  buffer4 << "Node Bits: " << bits;
  m.addDetail( buffer4.str() );
  return m;
}

class TriLagrangeShapeTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(TriLagrangeShapeTest);

    CPPUNIT_TEST(test_coeff_corners);
    CPPUNIT_TEST(test_coeff_edges);
    CPPUNIT_TEST(test_coeff_center);

    CPPUNIT_TEST(test_deriv_corners);
    CPPUNIT_TEST(test_deriv_edges);
    CPPUNIT_TEST(test_deriv_center);
    
    CPPUNIT_TEST_SUITE_END();
  
    TriLagrangeShape sf;
    
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



CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TriLagrangeShapeTest, "TriLagrangeShapeTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TriLagrangeShapeTest, "Unit");

double N0( double r, double   ) { return r*(2*r - 1); }
double N1( double  , double s ) { return s*(2*s - 1); }
double N2( double r, double s ) { double t = 1 - r - s; return t*(2*t - 1); }
double N3( double r, double s ) { return 4*r*s; }
double N4( double r, double s ) { double t = 1 - r - s; return 4*s*t; }
double N5( double r, double s ) { double t = 1 - r - s; return 4*r*t; }
double dN0dr( double r, double   ) { return 4*r - 1; }
double dN0ds( double  , double   ) { return 0; }
double dN1dr( double  , double   ) { return 0; }
double dN1ds( double  , double s ) { return 4*s - 1; }
double dN2dr( double r, double s ) { return 4*r + 4*s - 3; }
double dN2ds( double r, double s ) { return 4*r + 4*s - 3; }
double dN3dr( double  , double s ) { return 4*s; }
double dN3ds( double r, double   ) { return 4*r; }
double dN4dr( double  , double s ) { return -4*s; }
double dN4ds( double r, double s ) { return 4*(1-r-2*s); }
double dN5dr( double r, double s ) { return 4*(1-s-2*r); }
double dN5ds( double r, double   ) { return -4*r; }

typedef double (*N_t)(double,double);
const N_t N[] = { &N0, &N1, &N2, &N3, &N4, &N5 };
const N_t dNdr[] = { &dN0dr, &dN1dr, &dN2dr, &dN3dr, &dN4dr, &dN5dr };
const N_t dNds[] = { &dN0ds, &dN1ds, &dN2ds, &dN3ds, &dN4ds, &dN5ds };

const double rs_corner[][2] = { {1, 0}, {0, 1}, {0, 0}};
const double rs_edge[][2] = { {0.5, 0.5}, {0.0, 0.5}, {0.5, 0.0}};
const double rs_mid[2] = { 1.0/3.0, 1.0/3.0 };

static void get_coeff( NodeSet nodeset, const double* rs, double* coeffs )
{
  for (int i = 0; i < 6; ++i) 
    coeffs[i] = (*N[i])(rs[0], rs[1]);
  for (int i = 0; i < 3; ++i) 
    if (!nodeset.mid_edge_node(i)) {
      coeffs[i]      += 0.5 * coeffs[i+3];
      coeffs[(i+1)%3] += 0.5 * coeffs[i+3];
      coeffs[i+3] = 0;
    }
}

static void get_derivs( NodeSet nodeset, const double* rs, double* derivs )
{
  for (int i = 0; i < 6; ++i) {
    derivs[2*i  ] = (*dNdr[i])(rs[0], rs[1]);
    derivs[2*i+1] = (*dNds[i])(rs[0], rs[1]);
  }
  for (int i = 0; i < 3; ++i) 
    if (!nodeset.mid_edge_node(i)) {
      derivs[2*i]        += 0.5 * derivs[2*i+6];
      derivs[2*i+1]      += 0.5 * derivs[2*i+7];
      int j = (i+1)%3;
      derivs[2*j]        += 0.5 * derivs[2*i+6];
      derivs[2*j+1]      += 0.5 * derivs[2*i+7];
      derivs[2*i+6] = 0;
      derivs[2*i+7] = 0;
    }
}

static void check_valid_indices( const size_t* vtx_in, size_t num_vtx )
{
  CPPUNIT_ASSERT( num_vtx < 7 );
  CPPUNIT_ASSERT( num_vtx > 2 );
  size_t vertices[6];
  std::copy( vtx_in, vtx_in+num_vtx, vertices );
  std::sort( vertices, vertices + num_vtx );
  for (unsigned i = 1; i < num_vtx; ++i) {
    CPPUNIT_ASSERT( vertices[i] != vertices[i-1] );
    CPPUNIT_ASSERT( vertices[i] < 6 );
  }
}

static void check_no_zeros( const MsqVector<2>* derivs, size_t num_vtx )
{
  for (unsigned i = 0; i < num_vtx; ++i) {
    double dr = derivs[i][0]; 
    double ds = derivs[i][1]; 
    CPPUNIT_ASSERT( (fabs(dr) > 1e-6) || (fabs(ds) > 1e-6) );
  }
}

static void compare_coefficients( const double* coeffs,
                                  const size_t* indices,
                                  const double* expected_coeffs,
                                  size_t num_coeff,
                                  unsigned loc, NodeSet bits )
{
    // find the location in the returned list for each node
  size_t revidx[6];
  double test_vals[6];
  for (size_t i = 0; i < 6; ++i) {
    revidx[i] = std::find( indices, indices+num_coeff, i ) - indices;
    test_vals[i] = (revidx[i] == num_coeff) ? 0.0 : coeffs[revidx[i]];
  }

    // Check that index list doesn't contain any nodes not actually
    // present in the element.
  CPPUNIT_ASSERT( bits.mid_edge_node(0) || (revidx[3] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(1) || (revidx[4] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(2) || (revidx[5] == num_coeff) );
    
    // compare expected and actual coefficient values
  ASSERT_VALUES_EQUAL( expected_coeffs[0], test_vals[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[1], test_vals[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[2], test_vals[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[3], test_vals[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[4], test_vals[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[5], test_vals[5], loc, bits );
}

static void compare_derivatives( const size_t* vertices,
                                 size_t num_vtx,
                                 const MsqVector<2>* derivs,
                                 const double* expected_derivs,
                                 unsigned loc, NodeSet bits )
{
  check_valid_indices( vertices, num_vtx );
  check_no_zeros( derivs, num_vtx );
  double expanded_derivs[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (unsigned i = 0; i < num_vtx; ++i) {
    expanded_derivs[2*vertices[i]  ] = derivs[i][0];
    expanded_derivs[2*vertices[i]+1] = derivs[i][1];
  }
  
  ASSERT_VALUES_EQUAL( expected_derivs[ 0], expanded_derivs[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 1], expanded_derivs[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 2], expanded_derivs[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 3], expanded_derivs[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 4], expanded_derivs[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 5], expanded_derivs[5], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 6], expanded_derivs[ 6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 7], expanded_derivs[ 7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 8], expanded_derivs[ 8], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[ 9], expanded_derivs[ 9], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[10], expanded_derivs[10], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[11], expanded_derivs[11], loc, bits );
}

void TriLagrangeShapeTest::test_corner_coeff( int corner, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[6];
  get_coeff( nodebits, rs_corner[corner], expected );
  
  double coeff[27];
  size_t num_coeff = 17, indices[27];
  sf.coefficients( 0, corner, nodebits, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, expected, num_coeff, corner, nodebits );
}

void TriLagrangeShapeTest::test_edge_coeff( int edge, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[6];
  get_coeff( nodebits, rs_edge[edge], expected );
  
  double coeff[27];
  size_t num_coeff = 17, indices[27];
  sf.coefficients( 1, edge, nodebits, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, expected, num_coeff, edge+3, nodebits );
}

void TriLagrangeShapeTest::test_mid_coeff( NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[6];
  get_coeff( nodebits, rs_mid, expected );
  
  double coeff[27];
  size_t num_coeff = 17, indices[27];
  sf.coefficients( 2, 0, nodebits, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, expected, num_coeff, 6, nodebits );
}

void TriLagrangeShapeTest::test_corner_derivs( int corner, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[12];
  get_derivs( nodebits, rs_corner[corner], expected );
  
  size_t n = 19, vertices[100];
  MsqVector<2> derivs[100];
  sf.derivatives( 0, corner, nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, corner, nodebits );
}

void TriLagrangeShapeTest::test_edge_derivs( int edge, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[12];
  get_derivs( nodebits, rs_edge[edge], expected );
  
  size_t n = 19, vertices[100];
  MsqVector<2> derivs[100];
  sf.derivatives( 1, edge, nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, edge+3, nodebits );
}

void TriLagrangeShapeTest::test_mid_derivs( NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[12];
  get_derivs( nodebits, rs_mid, expected );
  
  size_t n = 19, vertices[100];
  MsqVector<2> derivs[100];
  sf.derivatives( 2, 0, nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, 6, nodebits );
}

void TriLagrangeShapeTest::test_coeff_corners()
{
  NodeSet ns;
  
  ns.clear();
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );

  ns.set_mid_edge_node(0);
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(1);
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );

  ns.set_mid_edge_node(0);
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(2);
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );

  ns.set_mid_edge_node(0);
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );

  ns.set_mid_edge_node(1);
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );

  ns.clear_mid_edge_node(0);
  test_corner_coeff( 0, ns );
  test_corner_coeff( 1, ns );
  test_corner_coeff( 2, ns );
}

void TriLagrangeShapeTest::test_coeff_edges()
{
  NodeSet ns;
  
  ns.clear();
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );

  ns.set_mid_edge_node(0);
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(1);
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );

  ns.set_mid_edge_node(0);
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(2);
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );

  ns.set_mid_edge_node(0);
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );

  ns.set_mid_edge_node(1);
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );

  ns.clear_mid_edge_node(0);
  test_edge_coeff( 0, ns );
  test_edge_coeff( 1, ns );
  test_edge_coeff( 2, ns );
}

void TriLagrangeShapeTest::test_coeff_center()
{
  NodeSet ns;
  
  ns.clear();
  test_mid_coeff( ns );
  ns.set_mid_edge_node(0);
  test_mid_coeff( ns );
  ns.clear();
  ns.set_mid_edge_node(1);
  test_mid_coeff( ns );
  ns.set_mid_edge_node(0);
  test_mid_coeff( ns );
  ns.clear();
  ns.set_mid_edge_node(2);
  test_mid_coeff( ns );
  ns.set_mid_edge_node(0);
  test_mid_coeff( ns );
  ns.set_mid_edge_node(1);
  test_mid_coeff( ns );
  ns.clear_mid_edge_node(0);
  test_mid_coeff( ns );
}

void TriLagrangeShapeTest::test_deriv_corners()
{
  NodeSet ns;
  
  ns.clear();
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );

  ns.set_mid_edge_node(0);
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(1);
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );

  ns.set_mid_edge_node(0);
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(2);
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );

  ns.set_mid_edge_node(0);
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );

  ns.set_mid_edge_node(1);
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );

  ns.clear_mid_edge_node(0);
  test_corner_derivs( 0, ns );
  test_corner_derivs( 1, ns );
  test_corner_derivs( 2, ns );
}

void TriLagrangeShapeTest::test_deriv_edges()
{
  NodeSet ns;
  
  ns.clear();
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );

  ns.set_mid_edge_node(0);
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(1);
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );

  ns.set_mid_edge_node(0);
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );

  ns.clear();
  ns.set_mid_edge_node(2);
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );

  ns.set_mid_edge_node(0);
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );

  ns.set_mid_edge_node(1);
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );

  ns.clear_mid_edge_node(0);
  test_edge_derivs( 0, ns );
  test_edge_derivs( 1, ns );
  test_edge_derivs( 2, ns );
}

void TriLagrangeShapeTest::test_deriv_center()
{
  NodeSet ns;
  
  ns.clear();
  test_mid_derivs( ns );
  ns.set_mid_edge_node(0);
  test_mid_derivs( ns );
  ns.clear();
  ns.set_mid_edge_node(1);
  test_mid_derivs( ns );
  ns.set_mid_edge_node(0);
  test_mid_derivs( ns );
  ns.clear();
  ns.set_mid_edge_node(2);
  test_mid_derivs( ns );
  ns.set_mid_edge_node(0);
  test_mid_derivs( ns );
  ns.set_mid_edge_node(1);
  test_mid_derivs( ns );
  ns.clear_mid_edge_node(0);
  test_mid_derivs( ns );
}
