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


/** \file LinearMappingFunctionTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MappingFunction.hpp"
#include "LinearHexahedron.hpp"
#include "LinearQuadrilateral.hpp"
#include "LinearTetrahedron.hpp"
#include "LinearTriangle.hpp"
#include "LinearPrism.hpp"
#include "LinearPyramid.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"

#include "UnitUtil.hpp"

#include <vector>
#include <algorithm>

using namespace Mesquite;
using namespace std;

class LinearMappingFunctionTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(LinearMappingFunctionTest);

    CPPUNIT_TEST(test_linear_hex_coeff_corners);
    CPPUNIT_TEST(test_linear_hex_coeff_edges);
    CPPUNIT_TEST(test_linear_hex_coeff_faces);
    CPPUNIT_TEST(test_linear_hex_coeff_center);

    CPPUNIT_TEST(test_linear_quad_coeff_corners);
    CPPUNIT_TEST(test_linear_quad_coeff_edges);
    CPPUNIT_TEST(test_linear_quad_coeff_faces);
    CPPUNIT_TEST(test_linear_quad_coeff_center);

    CPPUNIT_TEST(test_linear_tet_coeff_corners);
    CPPUNIT_TEST(test_linear_tet_coeff_edges);
    CPPUNIT_TEST(test_linear_tet_coeff_faces);
    CPPUNIT_TEST(test_linear_tet_coeff_center);

    CPPUNIT_TEST(test_linear_tri_coeff_corners);
    CPPUNIT_TEST(test_linear_tri_coeff_edges);
    CPPUNIT_TEST(test_linear_tri_coeff_faces);
    CPPUNIT_TEST(test_linear_tri_coeff_center);

    CPPUNIT_TEST(test_linear_prism_coeff_corners);
    CPPUNIT_TEST(test_linear_prism_coeff_edges);
    CPPUNIT_TEST(test_linear_prism_coeff_faces);
    CPPUNIT_TEST(test_linear_prism_coeff_center);

    CPPUNIT_TEST(test_linear_pyr_coeff_corners);
    CPPUNIT_TEST(test_linear_pyr_coeff_edges);
    CPPUNIT_TEST(test_linear_pyr_coeff_faces);
    CPPUNIT_TEST(test_linear_pyr_coeff_center);

    CPPUNIT_TEST(test_linear_hex_deriv_corners);
    CPPUNIT_TEST(test_linear_hex_deriv_edges);
    CPPUNIT_TEST(test_linear_hex_deriv_faces);
    CPPUNIT_TEST(test_linear_hex_deriv_center);

    CPPUNIT_TEST(test_linear_quad_deriv_corners);
    CPPUNIT_TEST(test_linear_quad_deriv_edges);
    CPPUNIT_TEST(test_linear_quad_deriv_faces);
    CPPUNIT_TEST(test_linear_quad_deriv_center);

    CPPUNIT_TEST(test_linear_tet_deriv_corners);
    CPPUNIT_TEST(test_linear_tet_deriv_edges);
    CPPUNIT_TEST(test_linear_tet_deriv_faces);
    CPPUNIT_TEST(test_linear_tet_deriv_center);

    CPPUNIT_TEST(test_linear_tri_deriv_corners);
    CPPUNIT_TEST(test_linear_tri_deriv_edges);
    CPPUNIT_TEST(test_linear_tri_deriv_faces);
    CPPUNIT_TEST(test_linear_tri_deriv_center);

    CPPUNIT_TEST(test_linear_prism_deriv_corners);
    CPPUNIT_TEST(test_linear_prism_deriv_edges);
    CPPUNIT_TEST(test_linear_prism_deriv_faces);
    CPPUNIT_TEST(test_linear_prism_deriv_center);

    CPPUNIT_TEST(test_linear_pyr_deriv_corners);
    CPPUNIT_TEST(test_linear_pyr_deriv_edges);
    CPPUNIT_TEST(test_linear_pyr_deriv_faces);
    CPPUNIT_TEST(test_linear_pyr_deriv_center);
    
    CPPUNIT_TEST_SUITE_END();
  
    LinearHexahedron hex;
    LinearQuadrilateral quad;
    LinearTetrahedron tet;
    LinearTriangle tri;
    LinearPrism prism;
    LinearPyramid pyr;
    
    static void hex_coeff( double xi[3], double coeff[8] );
    static void tet_coeff( double xi[3], double coeff[4] );
    static void quad_coeff( double xi[2], double coeff[4] );
    static void tri_coeff( double xi[2], double coeff[3] );
    static void prism_coeff( double xi[3], double coeff[3] );
    static void pyr_coeff( double xi[3], double coeff[3] );
    
    static void hex_deriv( double xi[3], double coeff_deriv[24] );
    static void tet_deriv( double xi[3], double coeff_deriv[12] );
    static void quad_deriv( double xi[2], double coeff_deriv[8] );
    static void tri_deriv( double xi[2], double coeff_deriv[6] );
    static void prism_deriv( double xi[2], double coeff_deriv[6] );
    static void pyr_deriv( double xi[2], double coeff_deriv[6] );
    
    typedef void (*map_func)(double*, double* );
    typedef void (MappingFunction::* mf_coeff)( unsigned, 
                                                unsigned, 
                                                double*,
                                                size_t&,
                                                MsqError& ) const;
    typedef void (MappingFunction::* mf_deriv)( unsigned, 
                                                unsigned, 
                                                size_t*,
                                                double*,
                                                size_t&,
                                                MsqError& ) const;
    
    void do_coeff_test( MappingFunction& mf, 
                        mf_coeff,
                        map_func mf2,
                        unsigned count,
                        double* xi );
    void do_deriv_test( MappingFunction& mf, 
                        mf_deriv,
                        map_func mf2,
                        unsigned count,
                        double* xi );
    void do_coeff_test_mid( MappingFunction& mf, 
                            map_func mf2,
                            double* xi );
    void do_deriv_test_mid( MappingFunction& mf, 
                            map_func mf2,
                             double* xi );
                  
    void do_test_fail( MappingFunction& mf, mf_coeff );
    void do_test_fail( MappingFunction& mf, mf_deriv );
    
    
    void xi_at_corners( EntityTopology type, double* xi, const int* corners );
    void xi_at_edges( EntityTopology type, double* xi, const int* corners );
    void xi_at_faces( EntityTopology type, double* xi, const int* corners );
    
  public:
  
    void setUp();
    void tearDown();
    

    void test_linear_hex_coeff_corners();
    void test_linear_hex_coeff_edges();
    void test_linear_hex_coeff_faces();
    void test_linear_hex_coeff_center();

    void test_linear_quad_coeff_corners();
    void test_linear_quad_coeff_edges();
    void test_linear_quad_coeff_faces();
    void test_linear_quad_coeff_center();

    void test_linear_tet_coeff_corners();
    void test_linear_tet_coeff_edges();
    void test_linear_tet_coeff_faces();
    void test_linear_tet_coeff_center();

    void test_linear_tri_coeff_corners();
    void test_linear_tri_coeff_edges();
    void test_linear_tri_coeff_faces();
    void test_linear_tri_coeff_center();

    void test_linear_prism_coeff_corners();
    void test_linear_prism_coeff_edges();
    void test_linear_prism_coeff_faces();
    void test_linear_prism_coeff_center();

    void test_linear_pyr_coeff_corners();
    void test_linear_pyr_coeff_edges();
    void test_linear_pyr_coeff_faces();
    void test_linear_pyr_coeff_center();

    void test_linear_hex_deriv_corners();
    void test_linear_hex_deriv_edges();
    void test_linear_hex_deriv_faces();
    void test_linear_hex_deriv_center();

    void test_linear_quad_deriv_corners();
    void test_linear_quad_deriv_edges();
    void test_linear_quad_deriv_faces();
    void test_linear_quad_deriv_center();

    void test_linear_tet_deriv_corners();
    void test_linear_tet_deriv_edges();
    void test_linear_tet_deriv_faces();
    void test_linear_tet_deriv_center();

    void test_linear_tri_deriv_corners();
    void test_linear_tri_deriv_edges();
    void test_linear_tri_deriv_faces();
    void test_linear_tri_deriv_center();

    void test_linear_prism_deriv_corners();
    void test_linear_prism_deriv_edges();
    void test_linear_prism_deriv_faces();
    void test_linear_prism_deriv_center();

    void test_linear_pyr_deriv_corners();
    void test_linear_pyr_deriv_edges();
    void test_linear_pyr_deriv_faces();
    void test_linear_pyr_deriv_center();
};



CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LinearMappingFunctionTest, "LinearMappingFunctionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LinearMappingFunctionTest, "Unit");

void LinearMappingFunctionTest::setUp() {}
void LinearMappingFunctionTest::tearDown() {}



/*******************************************************************************
 *                             Xi values at element corners
 *******************************************************************************/

const int HexSign[8][3] = { { -1, -1, -1 },
                            {  1, -1, -1 },
                            {  1,  1, -1 },
                            { -1,  1, -1 },
                            { -1, -1,  1 },
                            {  1, -1,  1 },
                            {  1,  1,  1 },
                            { -1,  1,  1 } };

const int QuadSign[4][2] = { { -1, -1 },
                             {  1, -1 },
                             {  1,  1 },
                             { -1,  1 } };


const int TetCorners[12] = { 0, 0, 0, 
                             1, 0, 0, 
                             0, 1, 0,
                             0, 0, 1 };
  
const int TriCorners[6] = { 0, 0,
                            1, 0, 
                            0, 1 };

const int PrismCorners[18] = { -1, 0, 0,
                               -1, 1, 0,
                               -1, 0, 1,
                                1, 0, 0,
                                1, 1, 0,
                                1, 0, 1 };
                                
/*******************************************************************************
 * Functions to calculate element xi at different locations, given the xi
 * at the corners.
 *******************************************************************************/

void LinearMappingFunctionTest::xi_at_corners( EntityTopology type, double* xi, const int* corners )
{
  unsigned d = TopologyInfo::dimension( type );
  unsigned c = TopologyInfo::corners( type );
  for (unsigned i = 0; i < c*d; ++i)
    xi[i] = corners[i];
}

void LinearMappingFunctionTest::xi_at_edges( EntityTopology type, double* xi, const int* corners )
{
  MsqError err;
  unsigned d = TopologyInfo::dimension( type );
  unsigned e = TopologyInfo::edges( type );
  for (unsigned i = 0; i < e; ++i)
  {
    const unsigned* vtx = TopologyInfo::edge_vertices( type, i, err );
    CPPUNIT_ASSERT(!err);
    for (unsigned j = 0; j < d; ++j)
      xi[d*i+j] = (corners[d*vtx[0]+j] + corners[d*vtx[1]+j])/2.0;
  }
}


void LinearMappingFunctionTest::xi_at_faces( EntityTopology type, double* xi, const int* corners )
{
  MsqError err;
  unsigned d = TopologyInfo::dimension( type );
  unsigned f = TopologyInfo::faces( type );
  for (unsigned i = 0; i < f; ++i)
  {
    unsigned c;
    const unsigned* vtx = TopologyInfo::face_vertices( type, i, c, err );
    CPPUNIT_ASSERT(!err);
    for (unsigned j = 0; j < d; ++j)
    {
      int sum = 0;
      for (unsigned k = 0; k < c; ++k)
        sum += corners[d*vtx[k]+j];
      xi[d*i+j] = (double)sum / c;
    }
  }
}




/*******************************************************************************
 *                 Test functions for mapping function coefficients
 *******************************************************************************/

void LinearMappingFunctionTest::test_linear_hex_coeff_corners()
{ 
  double xi[24];
  xi_at_corners( HEXAHEDRON, xi, &HexSign[0][0] );
  do_coeff_test( hex, &MappingFunction::coefficients_at_corner, hex_coeff, 8, xi );
}
  
void LinearMappingFunctionTest::test_linear_hex_coeff_edges()
{
  double xi[36];
  xi_at_edges( HEXAHEDRON, xi, &HexSign[0][0] );
  do_coeff_test( hex, &MappingFunction::coefficients_at_mid_edge, hex_coeff, 12, xi );
}

void LinearMappingFunctionTest::test_linear_hex_coeff_faces()
{
  double xi[18];
  xi_at_faces( HEXAHEDRON, xi, &HexSign[0][0] );
  do_coeff_test( hex, &MappingFunction::coefficients_at_mid_face, hex_coeff, 6, xi );
}

void LinearMappingFunctionTest::test_linear_hex_coeff_center()
{
  double xi[3] = { 0, 0, 0 };
  do_coeff_test_mid( hex, hex_coeff, xi );
}

void LinearMappingFunctionTest::test_linear_quad_coeff_corners()
{ 
  double xi[8];
  xi_at_corners( QUADRILATERAL, xi, &QuadSign[0][0] );
  do_coeff_test( quad, &MappingFunction::coefficients_at_corner, quad_coeff, 4, xi );
}

void LinearMappingFunctionTest::test_linear_quad_coeff_edges()
{
  double xi[8];
  xi_at_edges( QUADRILATERAL, xi, &QuadSign[0][0] );
  do_coeff_test( quad, &MappingFunction::coefficients_at_mid_edge, quad_coeff, 4, xi );
}

void LinearMappingFunctionTest::test_linear_quad_coeff_faces()
{
  do_test_fail( quad, &MappingFunction::coefficients_at_mid_face );
}

void LinearMappingFunctionTest::test_linear_quad_coeff_center()
{
  double xi[2] = { 0, 0 };
  do_coeff_test_mid( quad, quad_coeff, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_corners()
{
  double xi[12];
  xi_at_corners( TETRAHEDRON, xi, TetCorners );
  do_coeff_test( tet, &MappingFunction::coefficients_at_corner, tet_coeff, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_edges()
{
  double xi[18];
  xi_at_edges( TETRAHEDRON, xi, TetCorners );
  do_coeff_test( tet, &MappingFunction::coefficients_at_mid_edge, tet_coeff, 6, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_faces()
{
  double xi[12];
  xi_at_faces( TETRAHEDRON, xi, TetCorners );
  do_coeff_test( tet, &MappingFunction::coefficients_at_mid_face, tet_coeff, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_center()
{
  double xi[3] = { 0.25, 0.25, 0.25 };
  do_coeff_test_mid( tet, tet_coeff, xi );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_corners()
{
  double xi[12];
  xi_at_corners( TRIANGLE, xi, TriCorners );
  do_coeff_test( tri, &MappingFunction::coefficients_at_corner, tri_coeff, 3, xi );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_edges()
{
  double xi[18];
  xi_at_edges( TRIANGLE, xi, TriCorners );
  do_coeff_test( tri, &MappingFunction::coefficients_at_mid_edge, tri_coeff, 3, xi );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_faces()
{
  do_test_fail( tri, &MappingFunction::coefficients_at_mid_face );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_center()
{
  double xi[2] = { 1./3, 1./3 };
  do_coeff_test_mid( tri, tri_coeff, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_corners()
{
  double xi[18];
  xi_at_corners( PRISM, xi, PrismCorners );
  do_coeff_test( prism, &MappingFunction::coefficients_at_corner, prism_coeff, 6, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_edges()
{
  double xi[27];
  xi_at_edges( PRISM, xi, PrismCorners );
  do_coeff_test( prism, &MappingFunction::coefficients_at_mid_edge, prism_coeff, 9, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_faces()
{
  double xi[15];
  xi_at_faces( PRISM, xi, PrismCorners );
  do_coeff_test( prism, &MappingFunction::coefficients_at_mid_face, prism_coeff, 5, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_center()
{
  double xi[3] = { 0, 1./3, 1./3 };
  do_coeff_test_mid( prism, prism_coeff, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_corners()
{
  double xi[15] = { -1, -1, -1,
                     1, -1, -1,
                     1,  1, -1,
                    -1,  1, -1,
                     0,  0,  1 };
  do_coeff_test( pyr, &MappingFunction::coefficients_at_corner, pyr_coeff, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_edges()
{
  double xi[24] = { 0, -1, -1,
                    1,  0, -1,
                    0,  1, -1,
                   -1,  0, -1,
                   -1, -1,  0,
                    1, -1,  0,
                    1,  1,  0,
                   -1,  1,  0 };
  do_coeff_test( pyr, &MappingFunction::coefficients_at_mid_edge, pyr_coeff, 8, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_faces()
{
  double xi[15] = { 0, -1, -1./3,
                    1,  0, -1./3,
                    0,  1, -1./3,
                   -1,  0, -1./3,
                    0,  0, -1 };
  do_coeff_test( pyr, &MappingFunction::coefficients_at_mid_face, pyr_coeff, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_center()
{
  double xi[3] = { 0, 0, -0.5 };
  do_coeff_test_mid( pyr, pyr_coeff, xi );
}



/*******************************************************************************
 *                 Test functions for mapping function derivatives
 *******************************************************************************/

void LinearMappingFunctionTest::test_linear_hex_deriv_corners()
{ 
  double xi[24];
  xi_at_corners( HEXAHEDRON, xi, &HexSign[0][0] );
  do_deriv_test( hex, &MappingFunction::derivatives_at_corner, hex_deriv, 8, xi );
}
  
void LinearMappingFunctionTest::test_linear_hex_deriv_edges()
{
  double xi[36];
  xi_at_edges( HEXAHEDRON, xi, &HexSign[0][0] );
  do_deriv_test( hex, &MappingFunction::derivatives_at_mid_edge, hex_deriv, 12, xi );
}

void LinearMappingFunctionTest::test_linear_hex_deriv_faces()
{
  double xi[18];
  xi_at_faces( HEXAHEDRON, xi, &HexSign[0][0] );
  do_deriv_test( hex, &MappingFunction::derivatives_at_mid_face, hex_deriv, 6, xi );
}

void LinearMappingFunctionTest::test_linear_hex_deriv_center()
{
  double xi[3] = { 0, 0, 0 };
  do_deriv_test_mid( hex, hex_deriv, xi );
}

void LinearMappingFunctionTest::test_linear_quad_deriv_corners()
{ 
  double xi[8];
  xi_at_corners( QUADRILATERAL, xi, &QuadSign[0][0] );
  do_deriv_test( quad, &MappingFunction::derivatives_at_corner, quad_deriv, 4, xi );
}

void LinearMappingFunctionTest::test_linear_quad_deriv_edges()
{
  double xi[8];
  xi_at_edges( QUADRILATERAL, xi, &QuadSign[0][0] );
  do_deriv_test( quad, &MappingFunction::derivatives_at_mid_edge, quad_deriv, 4, xi );
}

void LinearMappingFunctionTest::test_linear_quad_deriv_faces()
{
  do_test_fail( quad, &MappingFunction::derivatives_at_mid_face );
}

void LinearMappingFunctionTest::test_linear_quad_deriv_center()
{
  double xi[2] = { 0, 0 };
  do_deriv_test_mid( quad, quad_deriv, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_corners()
{
  double xi[12];
  xi_at_corners( TETRAHEDRON, xi, TetCorners );
  do_deriv_test( tet, &MappingFunction::derivatives_at_corner, tet_deriv, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_edges()
{
  double xi[18];
  xi_at_edges( TETRAHEDRON, xi, TetCorners );
  do_deriv_test( tet, &MappingFunction::derivatives_at_mid_edge, tet_deriv, 6, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_faces()
{
  double xi[12];
  xi_at_faces( TETRAHEDRON, xi, TetCorners );
  do_deriv_test( tet, &MappingFunction::derivatives_at_mid_face, tet_deriv, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_center()
{
  double xi[3] = { 0.25, 0.25, 0.25 };
  do_deriv_test_mid( tet, tet_deriv, xi );
}
  
void LinearMappingFunctionTest::test_linear_tri_deriv_corners()
{
  double xi[12];
  xi_at_corners( TRIANGLE, xi, TriCorners );
  do_deriv_test( tri, &MappingFunction::derivatives_at_corner, tri_deriv, 3, xi );
}

void LinearMappingFunctionTest::test_linear_tri_deriv_edges()
{
  double xi[18];
  xi_at_edges( TRIANGLE, xi, TriCorners );
  do_deriv_test( tri, &MappingFunction::derivatives_at_mid_edge, tri_deriv, 3, xi );
}

void LinearMappingFunctionTest::test_linear_tri_deriv_faces()
{
  do_test_fail( tri, &MappingFunction::derivatives_at_mid_face );
}

void LinearMappingFunctionTest::test_linear_tri_deriv_center()
{
  double xi[2] = { 1./3, 1./3 };
  do_deriv_test_mid( tri, tri_deriv, xi );
}
   
void LinearMappingFunctionTest::test_linear_prism_deriv_corners()
{
  double xi[18];
  xi_at_corners( PRISM, xi, PrismCorners );
  do_deriv_test( prism, &MappingFunction::derivatives_at_corner, prism_deriv, 6, xi );
}

void LinearMappingFunctionTest::test_linear_prism_deriv_edges()
{
  double xi[27];
  xi_at_edges( PRISM, xi, PrismCorners );
  do_deriv_test( prism, &MappingFunction::derivatives_at_mid_edge, prism_deriv, 9, xi );
}

void LinearMappingFunctionTest::test_linear_prism_deriv_faces()
{
  double xi[15];
  xi_at_faces( PRISM, xi, PrismCorners );
  do_deriv_test( prism, &MappingFunction::derivatives_at_mid_face, prism_deriv, 5, xi );
}

void LinearMappingFunctionTest::test_linear_prism_deriv_center()
{
  double xi[3] = { 0, 1./3, 1./3 };
  do_deriv_test_mid( prism, prism_deriv, xi );
}
  


void LinearMappingFunctionTest::test_linear_pyr_deriv_corners()
{
  double xi[15] = { -1, -1, -1,
                     1, -1, -1,
                     1,  1, -1,
                    -1,  1, -1,
                     0,  0,  1 };
  do_deriv_test( pyr, &MappingFunction::derivatives_at_corner, pyr_deriv, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_deriv_edges()
{
  double xi[24] = { 0, -1, -1,
                    1,  0, -1,
                    0,  1, -1,
                   -1,  0, -1,
                   -1, -1,  0,
                    1, -1,  0,
                    1,  1,  0,
                   -1,  1,  0 };
  do_deriv_test( pyr, &MappingFunction::derivatives_at_mid_edge, pyr_deriv, 8, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_deriv_faces()
{
  double xi[15] = { 0, -1, -1./3,
                    1,  0, -1./3,
                    0,  1, -1./3,
                   -1,  0, -1./3,
                    0,  0, -1 };
  do_deriv_test( pyr, &MappingFunction::derivatives_at_mid_face, pyr_deriv, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_deriv_center()
{
  double xi[3] = { 0, 0, -0.5 };
  do_deriv_test_mid( pyr, pyr_deriv, xi );
}


/*******************************************************************************
 *        Mapping function implementations to compare with
 *******************************************************************************/
    
void LinearMappingFunctionTest::hex_coeff( double xi[3], double coeff[8] )
{
  for (unsigned i = 0; i < 8; ++i)
  {
    coeff[i] = 0.125;
    for (unsigned j = 0; j < 3; ++j)
      coeff[i] *= 1 + HexSign[i][j] * xi[j];
  }
}

void LinearMappingFunctionTest::tet_coeff( double xi[3], double coeff[4] )
{
  coeff[0] = 1 - xi[0] - xi[1] - xi[2];
  coeff[1] = xi[0];
  coeff[2] = xi[1];
  coeff[3] = xi[2];
}

void LinearMappingFunctionTest::quad_coeff( double xi[2], double coeff[4] )
{
  for (unsigned i = 0; i < 4; ++i)
  {
    coeff[i] = 0.25;
    for (unsigned j = 0; j < 2; ++j)
      coeff[i] *= 1 + QuadSign[i][j] * xi[j];
  }
}

void LinearMappingFunctionTest::tri_coeff( double xi[2], double coeff[3] )
{
  coeff[0] = 1 - xi[0] - xi[1];
  coeff[1] = xi[0];
  coeff[2] = xi[1];
}

void LinearMappingFunctionTest::prism_coeff( double xi[3], double coeff[6] )
{
  coeff[0] = 0.5 * (1 - xi[0]) * (1 - xi[1] - xi[2]);
  coeff[1] = 0.5 * (1 - xi[0]) * xi[1];
  coeff[2] = 0.5 * (1 - xi[0]) * xi[2];
  coeff[3] = 0.5 * (1 + xi[0]) * (1 - xi[1] - xi[2]);
  coeff[4] = 0.5 * (1 + xi[0]) * xi[1];
  coeff[5] = 0.5 * (1 + xi[0]) * xi[2];
}

void LinearMappingFunctionTest::pyr_coeff( double xi[3], double coeff[5] )
{
  coeff[0] = 0.125 * (1 - xi[0]) * (1 - xi[1]) * (1 - xi[2]);
  coeff[1] = 0.125 * (1 + xi[0]) * (1 - xi[1]) * (1 - xi[2]);
  coeff[2] = 0.125 * (1 + xi[0]) * (1 + xi[1]) * (1 - xi[2]);
  coeff[3] = 0.125 * (1 - xi[0]) * (1 + xi[1]) * (1 - xi[2]);
  coeff[4] = 0.500 *                             (1 + xi[2]);
}

/*******************************************************************************
 *        Mapping function derivatives to cmopare with
 *******************************************************************************/

void LinearMappingFunctionTest::hex_deriv( double xi[3], double coeff[24] )
{
  for (unsigned i = 0; i < 8; ++i)
    for (unsigned j = 0; j < 3; ++j)
      coeff[3*i+j] = 0.125 * HexSign[i][j] *
                     (1 + HexSign[i][(j+1)%3] * xi[(j+1)%3]) *
                     (1 + HexSign[i][(j+2)%3] * xi[(j+2)%3]);
}

void LinearMappingFunctionTest::tet_deriv( double*, double coeff[12] )
{
  static const double derivs[] = {-1,-1,-1,
                                   1, 0, 0, 
                                   0, 1, 0,
                                   0, 0, 1 };
  memcpy( coeff, derivs, sizeof(derivs) );
}

void LinearMappingFunctionTest::quad_deriv( double xi[2], double coeff[8] )
{
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 2; ++j)
      coeff[2*i+j] = 0.25 * QuadSign[i][j] * (1 + QuadSign[i][1-j] * xi[1-j]);
}

void LinearMappingFunctionTest::tri_deriv( double*, double coeff[6] )
{
  static const double derivs[] = {-1,-1,
                                   1, 0, 
                                   0, 1 };
  memcpy( coeff, derivs, sizeof(derivs) );
}

void LinearMappingFunctionTest::prism_deriv( double xi[3], double coeff[18] )
{
  coeff[ 0] = -0.5 * (1.0 - xi[1] - xi[2]);
  coeff[ 1] = -0.5 * (1.0 - xi[0]);
  coeff[ 2] = -0.5 * (1.0 - xi[0]);
  
  coeff[ 3] = -0.5 * xi[1];
  coeff[ 4] =  0.5 * (1.0 - xi[0]);
  coeff[ 5] =  0.0;
  
  coeff[ 6] = -0.5 * xi[2];
  coeff[ 7] =  0.0;
  coeff[ 8] =  0.5 * (1.0 - xi[0]);
  
  coeff[ 9] =  0.5 * (1.0 - xi[1] - xi[2]);
  coeff[10] = -0.5 * (1.0 + xi[0]);
  coeff[11] = -0.5 * (1.0 + xi[0]);
  
  coeff[12] =  0.5 * xi[1];
  coeff[13] =  0.5 * (1.0 + xi[0]);
  coeff[14] =  0.0;
  
  coeff[15] =  0.5 * xi[2];
  coeff[16] =  0.0;
  coeff[17] =  0.5 * (1.0 + xi[0]);
}

void LinearMappingFunctionTest::pyr_deriv( double xi[3], double coeff[15] )
{
  static int c[4][2] = { { -1, -1 },
                         {  1, -1 },
                         {  1,  1 },
                         { -1,  1 } };
  
  for (int n = 0; n < 4; ++n) {
    coeff[3*n  ] =  0.125 * c[n][0] * (1 + c[n][1] * xi[1]) * (1 - xi[2]);
    coeff[3*n+1] =  0.125 * c[n][1] * (1 + c[n][0] * xi[0]) * (1 - xi[2]);
    coeff[3*n+2] = -0.125 * (1 + c[n][0] * xi[0]) * (1 + c[n][1] * xi[1]);
  }
  coeff[12] = 0.0;
  coeff[13] = 0.0;
  coeff[14] = 0.5;
}

/*******************************************************************************
 *         Some utlity functions for creating CppUnit error messages
 *******************************************************************************/

static string itostr( int i )
{
  char buffer[32];
  sprintf(buffer, "%d", i );
  return buffer;
}

static string dtostr( double i )
{
  char buffer[32];
  sprintf(buffer, "%g", i );
  return buffer;
}



/*******************************************************************************
 *         Actual test imlplementation (common code for many tests)
 *******************************************************************************/
    
void LinearMappingFunctionTest::do_coeff_test( MappingFunction& mf, 
                                               mf_coeff func,
                                               map_func mf2,
                                               unsigned count,
                                               double* xi )
{
    // make sure it fails if passed a nonlinear element
  MsqError err;
  double coeff[100];
  size_t num_coeff = 100;
  (mf.*func)( 0, 1, coeff, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
  
    // get number of vertices in element
  const unsigned n = TopologyInfo::corners( mf.element_topology() );
  const unsigned d = TopologyInfo::dimension( mf.element_topology() );
  
    // compare coefficients at each location
  vector<double> comp(n);
  for (unsigned i = 0; i < count; ++i)
  {
    num_coeff = 101;
    (mf.*func)( i, 0, coeff, num_coeff, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(num_coeff, (size_t)n);
    
    mf2( xi, &comp[0] );
    string xi_str;
    for (unsigned j = 0; j < d; ++j) {
      xi_str += !j ? "(" : ", ";
      xi_str += dtostr(xi[j]);
    }
    xi_str += ")";
    xi += d;
    
    for (unsigned j = 0; j < n; ++j)
    {
      CppUnit::Message message( "Coefficients do not match." );
      message.addDetail( string("Entity:             ") + itostr( i ) );
      message.addDetail( string("Coefficient number: ") + itostr( j ) );
      message.addDetail( string("Xi:             ") + xi_str );
      message.addDetail( string("Expected value: ") + dtostr( comp[j] ) );
      message.addDetail( string("Actual value:   ") + dtostr( coeff[j] ) );
      ASSERT_MESSAGE( message, fabs(comp[j]-coeff[j]) < DBL_EPSILON );
    }
  }
}


void LinearMappingFunctionTest::do_deriv_test( MappingFunction& mf, 
                                               mf_deriv func,
                                               map_func mf2,
                                               unsigned count,
                                               double* xi )
{
    // make sure it fails if passed a nonlinear element
  MsqError err;
  double coeff[100];
  size_t verts[100], num_coeff = 37;
  (mf.*func)( 0, 1, verts, coeff, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
  
    // get number of vertices in element
  const unsigned n = TopologyInfo::corners( mf.element_topology() );
  const unsigned d = TopologyInfo::dimension( mf.element_topology() );
  
    // compare coefficients at each location
  vector<double> comp(n*d);
  for (unsigned i = 0; i < count; ++i)
  {
    num_coeff = 33;
    (mf.*func)( i, 0, verts, coeff, num_coeff, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( num_coeff > 0 );
    CPPUNIT_ASSERT( num_coeff <= n );
    
    mf2( xi, &comp[0] );
    string xi_str;
    for (unsigned j = 0; j < d; ++j) {
      xi_str += !j ? "(" : ", ";
      xi_str += dtostr(xi[j]);
    }
    xi_str += ")";
    xi += d;
    
    for (unsigned j = 0; j < num_coeff; ++j)
    {
      bool all_zero = true;
      for (unsigned k = 0; k < d; ++k)
      {
        CppUnit::Message message( "Coefficient derivatives do not match." );
        message.addDetail( string("Entity:             ") + itostr( i ) );
        message.addDetail( string("Coefficient number: ") + itostr( j ) );
        message.addDetail( string("Xi:             ") + xi_str );
        message.addDetail( string("Axis:           ") + itostr( k ) );
        message.addDetail( string("Expected value: ") + dtostr( comp[d*verts[j]+k] ) );
        message.addDetail( string("Actual value:   ") + dtostr( coeff[d*j+k] ) );
        ASSERT_MESSAGE( message, fabs(comp[d*verts[j]+k]-coeff[d*j+k]) < DBL_EPSILON );
        if (fabs(coeff[d*j+k]) > DBL_EPSILON)
          all_zero = false;
      }

        // if vertex has all zero values, it shouldn't have been in the
        // vertex list at all, as the Jacobian will not depend on that vertex.
      CPPUNIT_ASSERT( !all_zero );
    }
    
      // If any vertex is not in the list, then its values must be zero.
    sort( verts, verts + num_coeff );
    for (unsigned j = 0; j < num_coeff; ++j)
      if (!binary_search( verts, verts+num_coeff, j ))
        for (unsigned k = 0; k < d; ++k)
        {
          CppUnit::Message message( "Missing coefficient derivatives." );
          message.addDetail( string("Entity:              ") + itostr( i ) );
          message.addDetail( string("Coefficient number:  ") + itostr( j ) );
          message.addDetail( string("Axis:                ") + itostr( k ) );
          message.addDetail( string("Expected derivative: ") + dtostr( comp[d*j+k] ) );
          ASSERT_MESSAGE( message, fabs(comp[d*j+k]) < DBL_EPSILON );
        }
  }
}
                        
                        
void LinearMappingFunctionTest::do_coeff_test_mid( MappingFunction& mf, 
                                                   map_func mf2,
                                                   double* xi )
{
    // make sure it fails if passed a nonlinear element
  MsqError err;
  double coeff[100];
  size_t num_coeff = 38;
  mf.coefficients_at_mid_elem( 1, coeff, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
  
    // get number of vertices in element
  const unsigned n = TopologyInfo::corners( mf.element_topology() );
  
    // compare coefficients at each location
  vector<double> comp(n);
  num_coeff = 37;
  mf.coefficients_at_mid_elem( 0, coeff, num_coeff, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(num_coeff, (size_t)n);

  mf2( xi, &comp[0] );
  for (unsigned j = 0; j < n; ++j)
  {
    CppUnit::Message message( "Coefficients do not match." );
    message.addDetail( string("Coefficient number: ") + itostr( j ) );
    message.addDetail( string("Expected value: ") + dtostr( comp[j] ) );
    message.addDetail( string("Actual value:   ") + dtostr( coeff[j] ) );
    ASSERT_MESSAGE( message, fabs(comp[j]-coeff[j]) < DBL_EPSILON );
  }
}




void LinearMappingFunctionTest::do_deriv_test_mid( MappingFunction& mf, 
                                                   map_func mf2,
                                                   double* xi )
{
    // make sure it fails if passed a nonlinear element
  MsqError err;
  double coeff[100];
  size_t verts[100], num_vtx = 27;
  mf.derivatives_at_mid_elem( 1, verts, coeff, num_vtx, err );
  CPPUNIT_ASSERT(err);
  err.clear();
  
    // get number of vertices in element
  const unsigned n = TopologyInfo::corners( mf.element_topology() );
  const unsigned d = TopologyInfo::dimension( mf.element_topology() );
  
    // compare coefficients at each location
  vector<double> comp(n*d);
  num_vtx = 27;
  mf.derivatives_at_mid_elem( 0, verts, coeff, num_vtx, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( num_vtx > 0 );
  CPPUNIT_ASSERT( num_vtx <= n );

  mf2( xi, &comp[0] );
  for (unsigned j = 0; j < num_vtx; ++j)
  {
      // check that we got the expected values
    bool all_zero = true;
    for (unsigned k = 0; k < d; ++k)
    {
      CppUnit::Message message( "Coefficient derivatives do not match." );
      message.addDetail( string("Coefficient number: ") + itostr( j ) );
      message.addDetail( string("Axis:               ") + itostr( k ) );
      message.addDetail( string("Expected value: ") + dtostr( comp[d*verts[j]+k] ) );
      message.addDetail( string("Actual value:   ") + dtostr( coeff[d*j+k] ) );
      ASSERT_MESSAGE( message, fabs(comp[d*verts[j]+k]-coeff[d*j+k]) < DBL_EPSILON );
      if (fabs(coeff[d*j+k]) > DBL_EPSILON)
        all_zero = false;
    }

      // if vertex has all zero values, it shouldn't have been in the
      // vertex list at all, as the Jacobian will not depend on that vertex.
    CPPUNIT_ASSERT( !all_zero );
  }

    // If any vertex is not in the list, then its values must be zero.
  sort( verts, verts+num_vtx );
  for (unsigned j = 0; j < n; ++j)
    if (!binary_search( verts, verts+num_vtx, j ))
      for (unsigned k = 0; k < d; ++k)
      {
        CppUnit::Message message( "Missing coefficient derivatives." );
        message.addDetail( string("Coefficient number:  ") + itostr( j ) );
        message.addDetail( string("Axis:                ") + itostr( k ) );
        message.addDetail( string("Expected derivative: ") + dtostr( comp[d*j+k] ) );
        ASSERT_MESSAGE( message, fabs(comp[d*j+k]) < DBL_EPSILON );
      }
  
    // for linear elements, the derivative at the center should
    // depend on all the vertices
  CPPUNIT_ASSERT_EQUAL( (size_t)n, num_vtx );
}
                  
void LinearMappingFunctionTest::do_test_fail( MappingFunction& mf, mf_coeff func )
{
    // make sure it fails if called
  MsqError err;
  double coeff[100];
  size_t num_coeff;
  (mf.*func)( 0, 0, coeff, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
}  

void LinearMappingFunctionTest::do_test_fail( MappingFunction& mf, mf_deriv func )
{
    // make sure it fails if called
  MsqError err;
  double coeff[100];
  size_t verts[100], num_coeff;
  (mf.*func)( 0, 0, verts, coeff, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
}

