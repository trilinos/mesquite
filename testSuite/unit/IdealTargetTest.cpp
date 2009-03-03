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


/** \file IdealTargetTest.cpp
 *  \brief Test the IdealTargetCalculator
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdealTargetCalculator.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "PlanarDomain.hpp"
#include "MappingFunction.hpp"
#include "Settings.hpp"
#include "SamplePoints.hpp"
#include "IdealElements.hpp"
#include "ElemSampleQM.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "UnitUtil.hpp"

using namespace Mesquite;

class IdealTargetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(IdealTargetTest);
  CPPUNIT_TEST (test_tri_corner);
  CPPUNIT_TEST (test_tri_edge);
  CPPUNIT_TEST (test_tri_center);
  CPPUNIT_TEST (test_hex_corner);
  CPPUNIT_TEST (test_hex_edge);
  CPPUNIT_TEST (test_hex_face);
  CPPUNIT_TEST (test_hex_center);
  CPPUNIT_TEST (test_tri_orientation);
  CPPUNIT_TEST (test_quad_orientation);
  CPPUNIT_TEST (test_plane_neg_z);
  CPPUNIT_TEST_SUITE_END();

public:
  IdealTargetTest();
  void test_tri_corner();
  void test_tri_edge();
  void test_tri_center();
  void test_hex_corner();
  void test_hex_edge();
  void test_hex_face();
  void test_hex_center();
  void test_tri_orientation();
  void test_quad_orientation();
  void test_plane_neg_z();
  
private:

  void get_calc_target( bool rotate, EntityTopology type, 
                        unsigned dim, unsigned num,
                        MsqMatrix<3,3>&, MsqMatrix<3,2>& );
                        
  void get_ideal_target( EntityTopology type, 
                        unsigned dim, unsigned num,
                        MsqMatrix<3,3>&, MsqMatrix<3,2>& );
                        
  void do_test( EntityTopology type, unsigned dim, unsigned num );
  
  void compare_rotated( const MsqMatrix<3,2>& unrotated, 
                        const MsqMatrix<3,2>& rotated );
  
  const Vector3D planeNorm;
  Settings settings;
  PlanarDomain planeDomain;
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(IdealTargetTest, "IdealTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(IdealTargetTest, "Unit");

IdealTargetTest::IdealTargetTest()
  : planeNorm( 1, 0.5, -0.8 ),
    planeDomain( planeNorm, Vector3D( 0, 0, 0 ) )
  {}

void IdealTargetTest::test_tri_corner()
{
  do_test( TRIANGLE, 0, 0 );
  do_test( TRIANGLE, 0, 1 );
  do_test( TRIANGLE, 0, 2 );
}

void IdealTargetTest::test_tri_edge()
{
  do_test( TRIANGLE, 1, 0 );
  do_test( TRIANGLE, 1, 1 );
  do_test( TRIANGLE, 1, 2 );
}

void IdealTargetTest::test_tri_center()
{
  do_test( TRIANGLE, 2, 0 );
}

void IdealTargetTest::test_hex_corner()
{
  do_test( HEXAHEDRON, 0, 0 );
  do_test( HEXAHEDRON, 0, 1 );
  do_test( HEXAHEDRON, 0, 6 );
  do_test( HEXAHEDRON, 0, 7 );
}

void IdealTargetTest::test_hex_edge()
{
  do_test( HEXAHEDRON, 1, 1 );
  do_test( HEXAHEDRON, 1, 4 );
  do_test( HEXAHEDRON, 1, 11 );
}

void IdealTargetTest::test_hex_face()
{
  do_test( HEXAHEDRON, 2, 0 );
  do_test( HEXAHEDRON, 2, 1 );
  do_test( HEXAHEDRON, 2, 2 );
  do_test( HEXAHEDRON, 2, 3 );
  do_test( HEXAHEDRON, 2, 4 );
  do_test( HEXAHEDRON, 2, 5 );
}

void IdealTargetTest::test_hex_center()
{
  do_test( HEXAHEDRON, 3, 0 );
}

void IdealTargetTest::test_tri_orientation()
{
  MsqMatrix<3,2> W, Wr;
  MsqMatrix<3,3> junk;

  get_calc_target( false, TRIANGLE, 0, 0, junk, W );
  get_calc_target(  true, TRIANGLE, 0, 0, junk, Wr);
  compare_rotated( W, Wr );

  get_calc_target( false, TRIANGLE, 2, 0, junk, W );
  get_calc_target(  true, TRIANGLE, 2, 0, junk, Wr);
  compare_rotated( W, Wr );
}

void IdealTargetTest::test_quad_orientation()
{
  MsqMatrix<3,2> W, Wr;
  MsqMatrix<3,3> junk;

  get_calc_target( false, QUADRILATERAL, 0, 1, junk, W );
  get_calc_target(  true, QUADRILATERAL, 0, 1, junk, Wr);
  compare_rotated( W, Wr );

  get_calc_target( false, QUADRILATERAL, 2, 0, junk, W );
  get_calc_target(  true, QUADRILATERAL, 2, 0, junk, Wr);
  compare_rotated( W, Wr );
}

void IdealTargetTest::get_calc_target( bool rotate, EntityTopology type, 
                                       unsigned dim, unsigned num,
                                       MsqMatrix<3,3>& w3, MsqMatrix<3,2>& w2 )
{
  MsqPrintError err( msq_stdio::cout );
  const int elem_dim = TopologyInfo::dimension(type);
  
    // create a patch -- actual coords and such don't really matter
  std::vector<double> coords( 24, 0.0 );
  const size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  PatchData pd;
  pd.fill( 8, &coords[0], 1, type, conn, 0, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  pd.attach_settings( &settings );
  pd.set_domain( &planeDomain );
  
  SamplePoints pts( true, true, true, true );
  IdealTargetCalculator tc( rotate );
  unsigned sam = ElemSampleQM::sample( dim, num );
  if (elem_dim == 2)
    tc.get_2D_target( pd, 0, &pts, sam, w2, err );
  else
    tc.get_3D_target( pd, 0, &pts, sam, w3, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
}

void IdealTargetTest::get_ideal_target( EntityTopology type,
                                        unsigned dim, unsigned num,
                                        MsqMatrix<3,3>& w3, MsqMatrix<3,2>& w2 )
{
  MsqPrintError err( msq_stdio::cout );
  const unsigned elem_dim = TopologyInfo::dimension(type);
  
    // get the target matrix for an ideal element
  size_t indices[100];
  size_t num_vtx;
  const Vector3D* coords = unit_element( type );
  Vector3D c[3];
  if (elem_dim == 2) {
    MsqVector<2> derivs[100];
    const MappingFunction2D* func = settings.get_mapping_function_2D( type );
    func->derivatives( dim, num, 0, indices, derivs, num_vtx, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));

    for (size_t i = 0; i < num_vtx; ++i) 
      for (unsigned j = 0; j < 2; ++j)
        c[j] += derivs[i][j] * coords[indices[i]];

    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 2; ++j)
        w2(i,j) = c[j][i];
  }
  else {
    MsqVector<3> derivs[100];
    const MappingFunction3D* func = settings.get_mapping_function_3D( type );
    func->derivatives( dim, num, 0, indices, derivs, num_vtx, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));

    for (size_t i = 0; i < num_vtx; ++i) 
      for (unsigned j = 0; j < 3; ++j)
        c[j] += derivs[i][j] * coords[indices[i]];

    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        w3(i,j) = c[j][i];
  }
  
}

void IdealTargetTest::do_test( EntityTopology type, unsigned dim, unsigned num )
{
  MsqMatrix<3,3> w3_calc, w3_exp;
  MsqMatrix<3,2> w2_calc, w2_exp;
  get_calc_target( false, type, dim, num, w3_calc, w2_calc );
  get_ideal_target( type, dim, num, w3_exp, w2_exp );
  if (TopologyInfo::dimension(type) == 2)
    ASSERT_MATRICES_EQUAL( w2_exp, w2_calc, 1e-9 );
  else
    ASSERT_MATRICES_EQUAL( w3_exp, w3_calc, 1e-9 );
}

void IdealTargetTest::compare_rotated( const MsqMatrix<3,2>& unrotated, 
                                       const MsqMatrix<3,2>& rotated )
{
  MsqMatrix<3,1> u0, u1, r0, r1;
  u0 = unrotated.column(0);
  u1 = unrotated.column(1);
  r0 = rotated.column(0);
  r1 = rotated.column(1);
  MsqMatrix<3,1> un = u0 * u1, rn = r0 * r1;
  un /= length(un);
  rn /= length(rn);
  
    // unrotated element should be in XY plane
  const double z[] = { 0, 0, 1 };
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,1>(z)), un, 1e-6 );
  
    // rotated element should be in plane of domain
  MsqMatrix<3,1> pn( planeNorm.to_array() );
  pn /= length(pn);
  ASSERT_MATRICES_EQUAL( pn, rn, 1e-6 );
  
    // other properties should be unaffected by rotation
  CPPUNIT_ASSERT_DOUBLES_EQUAL( length(u0), length(r0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( length(u1), length(r1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( u0%u1, r0%r1, 1e-6 );
}

void IdealTargetTest::test_plane_neg_z()
{
  const Vector3D neg_z(0,0,-1);
  PlanarDomain dom(neg_z,Vector3D(0,0,0));
  MsqPrintError err( msq_stdio::cout );
  
    // create a patch -- actual coords and such don't really matter
  std::vector<double> coords( 24, 0.0 );
  const size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  PatchData pd;
  pd.fill( 8, &coords[0], 1, QUADRILATERAL, conn, 0, err );
  ASSERT_NO_ERROR(err);
  pd.attach_settings( &settings );
  pd.set_domain( &dom );
  
  SamplePoints pts( true, true, true, true );
  IdealTargetCalculator tc( true );
  
  unsigned sam = ElemSampleQM::sample( 0, 0) ;
  MsqMatrix<3,2> W;
  tc.get_2D_target( pd, 0, &pts, sam, W, err );
  ASSERT_NO_ERROR(err);
  
  Vector3D c1( W(0,0), W(1,0), W(2,0) );
  Vector3D c2( W(0,1), W(1,1), W(2,1) );
  Vector3D n = c1 * c2;
  n.normalize();
  CPPUNIT_ASSERT_VECTORS_EQUAL( neg_z, n, 1e-9 );
}

 
