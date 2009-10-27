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


/** \file LVQDTargetTest.cpp
 *  \brief unit tests for LVQDTargetCalculator class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LVQDTargetCalculator.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "UnitUtil.hpp"
#include "SizeFactor.hpp"
#include "OrientFactor.hpp"
#include "SkewFactor.hpp"
#include "AspectFactor.hpp"
#include "PatchDataInstances.hpp"

#include <iostream>

using namespace Mesquite;

class LVQDTargetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(LVQDTargetTest);

    // Test results of SkewFactor, OrientFactor, etc.
  CPPUNIT_TEST (test_factor_skew_2D);
  CPPUNIT_TEST (test_factor_skew_3D);
  CPPUNIT_TEST (test_factor_size_2D);
  CPPUNIT_TEST (test_factor_size_3D);
  CPPUNIT_TEST (test_factor_orient_2D);
  CPPUNIT_TEST (test_factor_orient_3D);
  CPPUNIT_TEST (test_factor_aspect_2D);
  CPPUNIT_TEST (test_factor_aspect_3D);

    // Test that LVQDTargetCalculator accepts a NULL component
    // and treats the corresponding value as I
  CPPUNIT_TEST (test_LVQD_default_is_I_2D);
  CPPUNIT_TEST (test_LVQD_default_is_I_3D);
    // Test that LVQDTargetCalculator returns the product of
    // its components.
  CPPUNIT_TEST (test_LVQD_product_2D);
  CPPUNIT_TEST (test_LVQD_product_3D);

    // Test that SkewFactor, OrientFactor, etc. correctly respond
    // to the underlying TargetCalculator flagging an error.
  CPPUNIT_TEST (test_factor_skew_error);
  CPPUNIT_TEST (test_factor_size_error);
  CPPUNIT_TEST (test_factor_orient_error);
  CPPUNIT_TEST (test_factor_aspect_error);
  CPPUNIT_TEST (test_LVQD_product_error);

    // Test that SkewFactor, OrientFactor, etc. correctly respond
    // to the underlying TargetCalculator returning 'invalid'.
  CPPUNIT_TEST (test_factor_skew_invalid);
  CPPUNIT_TEST (test_factor_size_invalid);
  CPPUNIT_TEST (test_factor_orient_invalid);
  CPPUNIT_TEST (test_factor_aspect_invalid);
  CPPUNIT_TEST (test_LVQD_product_invalid);

    // Test that SkewFactor, OrientFactor, etc. correctly respond
    // to the underlying TargetCalculator returning a zero-valued matrix.
  CPPUNIT_TEST (test_factor_skew_zero);
  CPPUNIT_TEST (test_factor_orient_zero);
  CPPUNIT_TEST (test_factor_aspect_zero);

  CPPUNIT_TEST_SUITE_END();

    // Define a few pre-factored matrices to use in tests.
    // Initialized by SetUp() method.
  MsqMatrix<3,3> V3D_Z45, V3D_X90, Q3D_45, D3D_123, I33;
  MsqMatrix<3,2> V2D_Z45, V2D_X90, I32;
  MsqMatrix<2,2> Q2D_45, D2D_21, I22;
    // PatchDatas for use in calls to 2D or 3D versions of functions.
  PatchData pd2D, pd3D;

    // Use LVQDTargetCalculator to calculate the product of the passed values.
  MsqMatrix<3,2> target( const double* L, 
                         const MsqMatrix<3,2>* V, 
                         const MsqMatrix<2,2>* Q, 
                         const MsqMatrix<2,2>* D );
    // Use LVQDTargetCalculator to calculate the product of the passed values.
  MsqMatrix<3,3> target( const double* L, 
                         const MsqMatrix<3,3>* V, 
                         const MsqMatrix<3,3>* Q, 
                         const MsqMatrix<3,3>* D );

    // Use SizeFactor to extract size component of M
  double L_factor( MsqMatrix<3,2> M );
  double L_factor( MsqMatrix<2,2> M );
  double L_factor( MsqMatrix<3,3> M );

    // Use OrientFactor to extract orientation component of M
  MsqMatrix<3,2> V_factor( MsqMatrix<3,2> M );
  MsqMatrix<3,3> V_factor( MsqMatrix<3,3> M );

    // Use SkewFactor to extract skew component of M
  MsqMatrix<2,2> Q_factor( MsqMatrix<3,2> M );
  MsqMatrix<2,2> Q_factor( MsqMatrix<2,2> M );
  MsqMatrix<3,3> Q_factor( MsqMatrix<3,3> M );

    // Use AspectFactor to extract aspect ratio component of M
  MsqMatrix<2,2> D_factor( MsqMatrix<3,2> M );
  MsqMatrix<2,2> D_factor( MsqMatrix<2,2> M );
  MsqMatrix<3,3> D_factor( MsqMatrix<3,3> M );
  
    // test that two matrices are a rotation of each other.
  static void check_is_rotation( MsqMatrix<2,2> A, MsqMatrix<2,2> B );
  static void check_is_rotation( MsqMatrix<3,3> A, MsqMatrix<3,3> B );
  
    // Test that SkewFactor, etc return 'invalid' and optionally flag an error.
  void check_skew_invalid( TargetCalculator& Wcalc, bool expect_error );
  void check_size_invalid( TargetCalculator& Wcalc, bool expect_error );
  void check_orient_invalid( TargetCalculator& Wcalc, bool expect_error );
  void check_aspect_invalid( TargetCalculator& Wcalc, bool expect_error );

public:

  void setUp();

    // Test results of SkewFactor, OrientFactor, etc.
  void test_factor_skew_2D();
  void test_factor_skew_3D();
  void test_factor_size_2D();
  void test_factor_size_3D();
  void test_factor_orient_2D();
  void test_factor_orient_3D();
  void test_factor_aspect_2D();
  void test_factor_aspect_3D();

    // Test that LVQDTargetCalculator accepts a NULL component
    // and treats the corresponding value as I
  void test_LVQD_default_is_I_2D();
  void test_LVQD_default_is_I_3D();
    // Test that LVQDTargetCalculator returns the product of
    // its components.
  void test_LVQD_product_2D();
  void test_LVQD_product_3D();

    // Test that SkewFactor, OrientFactor, etc. correctly respond
    // to the underlying TargetCalculator flagging an error.
  void test_factor_skew_error();
  void test_factor_size_error();
  void test_factor_orient_error();
  void test_factor_aspect_error();
  void test_LVQD_product_error();

    // Test that SkewFactor, OrientFactor, etc. correctly respond
    // to the underlying TargetCalculator returning 'invalid'.
  void test_factor_skew_invalid();
  void test_factor_size_invalid();
  void test_factor_orient_invalid();
  void test_factor_aspect_invalid();
  void test_LVQD_product_invalid();

    // Test that SkewFactor, OrientFactor, etc. correctly respond
    // to the underlying TargetCalculator returning a zero-valued matrix.
  void test_factor_skew_zero();
  void test_factor_orient_zero();
  void test_factor_aspect_zero();

    // Helper class: return constant values for target matrices.
  class ConstantTarget : public TargetCalculator
  {
    private:
      MsqMatrix<3,3> target3D;
      MsqMatrix<3,2> target2D;
      bool have3D, have2D;
      bool flagError;

    public:
      ConstantTarget( MsqMatrix<3,3> val3D, MsqMatrix<3,2> val2D ) 
        : target3D(val3D), target2D(val2D), have3D(true), have2D(true) 
        {}
      ConstantTarget( MsqMatrix<3,3> val3D, MsqMatrix<2,2> val2D ) 
        : target3D(val3D), target2D(0.0), have3D(true), have2D(true) 
        { target2D.set_row(0,val2D.row(0)); target2D.set_row(1,val2D.row(1)); }
      ConstantTarget( double C )
        : target3D(C), target2D(C), have3D(true), have2D(true) 
        {}
      ConstantTarget( MsqMatrix<3,3> val3D ) 
        : target3D(val3D), target2D(0.0), have3D(true), have2D(false) 
        {}
      ConstantTarget( MsqMatrix<3,2> val2D ) 
        : target3D(0.0), target2D(val2D), have3D(false), have2D(true) 
        {}
      ConstantTarget( MsqMatrix<2,2> val2D ) 
        : target3D(0.0), target2D(0.0), have3D(false), have2D(true) 
        { target2D.set_row(0,val2D.row(0)); target2D.set_row(1,val2D.row(1)); }

      virtual bool get_3D_target( PatchData&, size_t, Sample, MsqMatrix<3,3>& result, MsqError& err )
        { CPPUNIT_ASSERT(have3D); result = target3D; return true; }

      virtual bool get_2D_target( PatchData&, size_t, Sample, MsqMatrix<3,2>& result, MsqError& err )
        { CPPUNIT_ASSERT(have2D); result = target2D; return true; }
  };

    // Helper class: return 'invalid' for target matrices, and
    // optionally flag an error.
  class TargetError : public TargetCalculator
  {
      bool flagError;

    public:
      TargetError( bool flag_error = true ) : flagError(flag_error) {}

      bool get_3D_target( PatchData&, size_t, Sample, MsqMatrix<3,3>&, MsqError& err)
        { if (flagError) MSQ_SETERR(err)(MsqError::INVALID_MESH); return false; }

      bool get_2D_target( PatchData&, size_t, Sample, MsqMatrix<3,2>&, MsqError& err )
        { if (flagError) MSQ_SETERR(err)(MsqError::INVALID_MESH); return false; }
  };
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LVQDTargetTest, "LVQDTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LVQDTargetTest, "Unit");


void LVQDTargetTest::setUp()
{
  const double cos45 = MSQ_SQRT_TWO/2.0;
  const double rotation_3D_Z45[9] = { cos45, -cos45, 0,
                                      cos45,  cos45, 0,
                                          0,      0, 1 };
  const double rotation_2D_Z45[6] = { cos45, -cos45,
                                      cos45,  cos45,
                                          0,      0 };
                                          
  const double rotation_3D_X90[9] = { 1,  0,  0,
                                      0,  0, -1,
                                      0,  1,  0 };
  const double rotation_2D_X90[6] = { 1,  0,
                                      0,  0,
                                      0,  1 };
                                      
  const double rc45 = sqrt(cos45);
  const double skew_2D_45[4] = { 1/rc45, rc45,
                                 0,      rc45 };
  const double skew_3D_45[9] = { 1, cos45, cos45,
                                 0, cos45, 1 - cos45,
                                 0,     0, sqrt(MSQ_SQRT_TWO-1) };
                                 
  const double aspect_2D_2x[4] = { MSQ_SQRT_TWO, 0,
                                   0,            MSQ_SQRT_TWO/2 };
  const double r6 = Mesquite::cbrt(1.0/6.0);
  const double aspect_3D_123[9] = { r6,    0,    0,
                                     0, 2*r6,    0,
                                     0,    0, 3*r6 };
                                     
  V3D_Z45 = MsqMatrix<3,3>(rotation_3D_Z45);
  V3D_X90 = MsqMatrix<3,3>(rotation_3D_X90);
  Q3D_45  = MsqMatrix<3,3>(skew_3D_45);
  Q3D_45  *= 1/Mesquite::cbrt(det(Q3D_45));
  D3D_123 = MsqMatrix<3,3>(aspect_3D_123);
  
  V2D_Z45 = MsqMatrix<3,2>(rotation_2D_Z45);
  V2D_X90 = MsqMatrix<3,2>(rotation_2D_X90);
  Q2D_45  = MsqMatrix<2,2>(skew_2D_45);
  D2D_21  = MsqMatrix<2,2>(aspect_2D_2x);
  
  I33 = MsqMatrix<3,3>(1.0);
  I32 = MsqMatrix<3,2>(1.0);
  I22 = MsqMatrix<2,2>(1.0);
  
  MsqError err;
  create_one_tet_patch( pd3D, err );
  ASSERT_NO_ERROR(err);
  create_one_tri_patch( pd2D, err );
  ASSERT_NO_ERROR(err);
};

void LVQDTargetTest::check_is_rotation( MsqMatrix<2,2> A, MsqMatrix<2,2> B )
{
  MsqMatrix<2,2> R = B * inverse(A);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( R(0,0), R(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( R(0,1), -R(1,0), 1e-6 );
}

void LVQDTargetTest::check_is_rotation( MsqMatrix<3,3> A, MsqMatrix<3,3> B )
{
  MsqMatrix<3,3> R = B * inverse(A);
  MsqMatrix<3,1> C0 = R.column(0);
  MsqMatrix<3,1> C1 = R.column(1);
  MsqMatrix<3,1> C2 = R.column(2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, length(C0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, length(C1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, length(C2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, C0 % C1, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, C0 % C2, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, C1 % C2, 1e-6 );
}


MsqMatrix<2,2> LVQDTargetTest::Q_factor( MsqMatrix<3,2> M )
{
  MsqError err;
  MsqMatrix<2,2> Q;
  ConstantTarget W( M );
  SkewFactor F( &W );
  bool v = F.get_skew_2D( pd2D, 0, Sample(0,0), Q, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return Q;
}
MsqMatrix<2,2> LVQDTargetTest::Q_factor( MsqMatrix<2,2> M )
{
  MsqError err;
  MsqMatrix<2,2> Q;
  ConstantTarget W( M );
  SkewFactor F( &W );
  bool v = F.get_skew_2D( pd2D, 0, Sample(0,0), Q, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return Q;
}
MsqMatrix<3,3> LVQDTargetTest::Q_factor( MsqMatrix<3,3> M )
{
  MsqError err;
  MsqMatrix<3,3> Q;
  ConstantTarget W( M );
  SkewFactor F( &W );
  bool v = F.get_skew_3D( pd3D, 0, Sample(0,0), Q, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return Q;
}

void LVQDTargetTest::test_factor_skew_2D()
{
    // test that, given only the skew portion of the
    // matrix, we get back the input
  check_is_rotation( I22,    Q_factor(I32) );
  check_is_rotation( Q2D_45, Q_factor(Q2D_45) );
  
    // test that scaling is removed
  check_is_rotation( I22, Q_factor(MsqMatrix<3,2>(5.0)) );
  check_is_rotation( Q2D_45, Q_factor(3*Q2D_45) );
  
    // test that rotation and aspect ratio changes are removed
  check_is_rotation( I22,    Q_factor(0.5*V2D_Z45*       D2D_21) );
  check_is_rotation( Q2D_45, Q_factor(1.5*V2D_Z45*Q2D_45*D2D_21) );
}  

void LVQDTargetTest::test_factor_skew_3D()
{
    // test that, given only the skew portion of the
    // matrix, we get back the input
  check_is_rotation( I33,    Q_factor(I33) );
  check_is_rotation( Q3D_45, Q_factor(Q3D_45) );
  
    // test that scaling is removed
  check_is_rotation( I33, Q_factor(MsqMatrix<3,3>(5.0)) );
  check_is_rotation( Q3D_45, Q_factor(3*Q3D_45) );
  
    // test that rotation and aspect ratio changes are removed
  check_is_rotation( I33,    Q_factor(0.5*V3D_Z45*       D3D_123) );
  check_is_rotation( Q3D_45, Q_factor(1.5*V3D_Z45*Q3D_45*D3D_123) );
}


double LVQDTargetTest::L_factor( MsqMatrix<3,2> M )
{
  MsqError err;
  ConstantTarget W( M );
  SizeFactor F( &W );
  double L;
  bool v = F.get_size( pd2D, 0, Sample(0,0), L, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return L;
}
double LVQDTargetTest::L_factor( MsqMatrix<2,2> M )
{
  MsqError err;
  ConstantTarget W( M );
  SizeFactor F( &W );
  double L;
  bool v = F.get_size( pd2D, 0, Sample(0,0), L, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return L;
}
double LVQDTargetTest::L_factor( MsqMatrix<3,3> M )
{
  MsqError err;
  ConstantTarget W( M );
  SizeFactor F( &W );
  double L;
  bool v = F.get_size( pd3D, 0, Sample(0,0), L, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return L;
}

void LVQDTargetTest::test_factor_size_2D()
{
    // check that we get 1.0 for all matrices with unit size
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, L_factor(I32), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, L_factor(V2D_Z45), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, L_factor(Q2D_45), 1e-8 );

    // scale matrices and test for expected scaled size
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, L_factor(MsqMatrix<3,2>(5.0)), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.3, L_factor(0.3*V2D_Z45), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.1, L_factor(1.1*Q2D_45), 1e-8 );
}

void LVQDTargetTest::test_factor_size_3D()
{
    // check that we get 1.0 for all matrices with unit size
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, L_factor(I32), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, L_factor(V3D_Z45), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, L_factor(Q3D_45), 1e-8 );

    // scale matrices and test for expected scaled size
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, L_factor(MsqMatrix<3,3>(5.0)), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.3, L_factor(0.3*V3D_Z45), 1e-8 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.1, L_factor(1.1*Q3D_45), 1e-8 );
}


MsqMatrix<3,2> LVQDTargetTest::V_factor( MsqMatrix<3,2> M )
{
  MsqError err;
  MsqMatrix<3,2> V;
  ConstantTarget W( M );
  OrientFactor F( &W );
  bool v = F.get_orient_2D( pd2D, 0, Sample(0,0), V, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return V;
}
MsqMatrix<3,3> LVQDTargetTest::V_factor( MsqMatrix<3,3> M )
{
  MsqError err;
  MsqMatrix<3,3> V;
  ConstantTarget W( M );
  OrientFactor F( &W );
  bool v = F.get_orient_3D( pd3D, 0, Sample(0,0), V, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return V;
}

void LVQDTargetTest::test_factor_orient_2D()
{
    // test that, given only the orientation portion of the
    // matrix, we get back the input
  ASSERT_MATRICES_EQUAL( I32,     V_factor(I32),     1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45, V_factor(V2D_Z45), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_X90, V_factor(V2D_X90), 1e-8 );
  
    // test that scaling is removed
  ASSERT_MATRICES_EQUAL( I32,     V_factor(5.0*I32),     1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45, V_factor(2.3*V2D_Z45), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_X90, V_factor(0.6*V2D_X90), 1e-8 );
  
    // test factorization of matrix containing all factors
  ASSERT_MATRICES_EQUAL( I32,     V_factor(0.3*I32    *Q2D_45*D2D_21), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45, V_factor(1.2*V2D_Z45*Q2D_45*D2D_21), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_X90, V_factor(2.7*V2D_X90*Q2D_45*D2D_21), 1e-8 );
}

void LVQDTargetTest::test_factor_orient_3D()
{
    // test that, given only the orientation portion of the
    // matrix, we get back the input
  ASSERT_MATRICES_EQUAL( I33,     V_factor(I33),     1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45, V_factor(V3D_Z45), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_X90, V_factor(V3D_X90), 1e-8 );
  
    // test that scaling is removed
  ASSERT_MATRICES_EQUAL( I33,     V_factor(5.0*I33),     1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45, V_factor(2.3*V3D_Z45), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_X90, V_factor(0.6*V3D_X90), 1e-8 );
  
    // test factorization of matrix containing all factors
  ASSERT_MATRICES_EQUAL( I33,     V_factor(0.3*        Q3D_45*D3D_123), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45, V_factor(1.2*V3D_Z45*Q3D_45*D3D_123), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_X90, V_factor(2.7*V3D_X90*Q3D_45*D3D_123), 1e-8 );
}


MsqMatrix<2,2> LVQDTargetTest::D_factor( MsqMatrix<3,2> M )
{
  MsqError err;
  MsqMatrix<2,2> D;
  ConstantTarget W( M );
  AspectFactor F( &W );
  bool v = F.get_aspect_2D( pd2D, 0, Sample(0,0), D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return D;
}
MsqMatrix<2,2> LVQDTargetTest::D_factor( MsqMatrix<2,2> M )
{
  MsqError err;
  MsqMatrix<2,2> D;
  ConstantTarget W( M );
  AspectFactor F( &W );
  bool v = F.get_aspect_2D( pd2D, 0, Sample(0,0), D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return D;
}
MsqMatrix<3,3> LVQDTargetTest::D_factor( MsqMatrix<3,3> M )
{
  MsqError err;
  MsqMatrix<3,3> D;
  ConstantTarget W( M );
  AspectFactor F( &W );
  bool v = F.get_aspect_3D( pd3D, 0, Sample(0,0), D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return D;
}

void LVQDTargetTest::test_factor_aspect_2D()
{
    // test that, given only the apect portion of the
    // matrix, we get back the input
  ASSERT_MATRICES_EQUAL( I22,    D_factor(I32),    1e-8 );
  ASSERT_MATRICES_EQUAL( D2D_21, D_factor(D2D_21), 1e-8 );
  
    // test that scaling is removed
  ASSERT_MATRICES_EQUAL( I22,    D_factor(MsqMatrix<3,2>(5.0)), 1e-8 );
  ASSERT_MATRICES_EQUAL( D2D_21, D_factor(3*D2D_21),            1e-8 );
  
    // test with composition of all factors
  ASSERT_MATRICES_EQUAL( I22,    D_factor(0.5*V2D_Z45*Q2D_45       ), 1e-8 );
  ASSERT_MATRICES_EQUAL( D2D_21, D_factor(1.5*V2D_Z45*Q2D_45*D2D_21), 1e-8 );
}

void LVQDTargetTest::test_factor_aspect_3D()

{
    // test that, given only the apect portion of the
    // matrix, we get back the input
  ASSERT_MATRICES_EQUAL( I33,     D_factor(I33),    1e-8 );
  ASSERT_MATRICES_EQUAL( D3D_123, D_factor(D3D_123), 1e-8 );
  
    // test that scaling is removed
  ASSERT_MATRICES_EQUAL( I33,     D_factor(5*I33),    1e-8 );
  ASSERT_MATRICES_EQUAL( D3D_123, D_factor(3*D3D_123), 1e-8 );
  
    // test with composition of all factors
  ASSERT_MATRICES_EQUAL( I33,     D_factor(0.5*V3D_Z45*Q3D_45        ), 1e-8 );
  ASSERT_MATRICES_EQUAL( D3D_123, D_factor(1.5*V3D_X90*Q3D_45*D3D_123), 1e-8 );
}

MsqMatrix<3,2> LVQDTargetTest::target( const double* L, 
                                       const MsqMatrix<3,2>* V, 
                                       const MsqMatrix<2,2>* Q, 
                                       const MsqMatrix<2,2>* D )
{
  ConstantTarget W_size  ( L ? *L : 1.0 );
  ConstantTarget W_orient( V ? *V : I32 );
  ConstantTarget W_skew  ( Q ? *Q : I22 );
  ConstantTarget W_aspect( D ? *D : I22 );
  SizeFactor   F_size  ( &W_size   );
  OrientFactor F_orient( &W_orient );
  SkewFactor   F_skew  ( &W_skew   );
  AspectFactor F_aspect( &W_aspect );
  LVQDTargetCalculator LVQD( L ? &F_size   : NULL,
                             V ? &F_orient : NULL,
                             Q ? &F_skew   : NULL,
                             D ? &F_aspect : NULL );
  MsqError err;
  MsqMatrix<3,2> W;
  bool v = LVQD.get_2D_target( pd2D, 0, Sample(0,0), W, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return W;
}
MsqMatrix<3,3> LVQDTargetTest::target( const double* L, 
                                       const MsqMatrix<3,3>* V, 
                                       const MsqMatrix<3,3>* Q, 
                                       const MsqMatrix<3,3>* D )
{
  ConstantTarget W_size  ( L ? *L : 1.0 );
  ConstantTarget W_orient( V ? *V : I33 );
  ConstantTarget W_skew  ( Q ? *Q : I33 );
  ConstantTarget W_aspect( D ? *D : I33 );
  SizeFactor   F_size  ( &W_size   );
  OrientFactor F_orient( &W_orient );
  SkewFactor   F_skew  ( &W_skew   );
  AspectFactor F_aspect( &W_aspect );
  LVQDTargetCalculator LVQD( L ? &F_size   : NULL,
                             V ? &F_orient : NULL,
                             Q ? &F_skew   : NULL,
                             D ? &F_aspect : NULL );
  MsqError err;
  MsqMatrix<3,3> W;
  bool v = LVQD.get_3D_target( pd3D, 0, Sample(0,0), W, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return W;
}

void LVQDTargetTest::test_LVQD_default_is_I_2D()
{
  double s = 3.2;
  MsqMatrix<3,2>* null_V = 0;
  ASSERT_MATRICES_EQUAL( I32, target(0,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,2>(s)), target(&s,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45, target(0,&V2D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V2D_Z45, target(&s,&V2D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*Q2D_45, target(0,&V2D_Z45,&Q2D_45,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*D2D_21, target(0,&V2D_Z45,0,&D2D_21), 1e-8 );
}

void LVQDTargetTest::test_LVQD_default_is_I_3D()
{
  double s = 2.6;
  MsqMatrix<3,3>* null_V = 0;
  ASSERT_MATRICES_EQUAL( I33, target(0,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(s)), target(&s,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45, target(0,&V3D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V3D_Z45, target(&s,&V3D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*Q3D_45,  target(0,&V3D_Z45,&Q3D_45, 0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*D3D_123, target(0,&V3D_Z45,0,&D3D_123), 1e-8 );
}

void LVQDTargetTest::test_LVQD_product_2D()
{
  double s = 3.2;
  double o = 1.0;
  ASSERT_MATRICES_EQUAL( I32, target(&o,&I32,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,2>(s)), target(&s,&I32,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45, target(&o,&V2D_Z45,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V2D_Z45, target(&s,&V2D_Z45,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*Q2D_45, target(&o,&V2D_Z45,&Q2D_45,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*D2D_21, target(&o,&V2D_Z45,&I22,&D2D_21), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V2D_Z45*Q2D_45*D2D_21, target(&s,&V2D_Z45,&Q2D_45,&D2D_21), 1e-8 );
}

void LVQDTargetTest::test_LVQD_product_3D()
{
  double s = 2.6;
  double o = 1.0;
  ASSERT_MATRICES_EQUAL( I33, target(&o,&I33,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(s)), target(&s,&I33,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45, target(&o,&V3D_Z45,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V3D_Z45, target(&s,&V3D_Z45,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*Q3D_45, target(&o,&V3D_Z45,&Q3D_45,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*D3D_123, target(&o,&V3D_Z45,&I33,&D3D_123), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V3D_Z45*Q3D_45*D3D_123, target(&s,&V3D_Z45,&Q3D_45,&D3D_123), 1e-8 );
}

  
void LVQDTargetTest::check_skew_invalid( TargetCalculator& Wcalc, 
                                         bool expect_error )
{
  SkewFactor Qcalc(&Wcalc);
  MsqError err;

  MsqMatrix<2,2> Q2D;
  bool valid = Qcalc.get_skew_2D( pd2D, 0, Sample(0,0), Q2D, err );
  CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
  CPPUNIT_ASSERT(!valid);

  MsqMatrix<3,3> Q3D;
  valid = Qcalc.get_skew_3D( pd3D, 0, Sample(0,0), Q3D, err );
  CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
  CPPUNIT_ASSERT(!valid);
}

void LVQDTargetTest::check_size_invalid( TargetCalculator& Wcalc, 
                                         bool expect_error )
{
  SizeFactor Lcalc(&Wcalc);
  MsqError err;

  double lambda;
  bool valid = Lcalc.get_size( pd2D, 0, Sample(0,0), lambda, err );
  CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
  CPPUNIT_ASSERT(!valid);

  valid = Lcalc.get_size( pd3D, 0, Sample(0,0), lambda, err );
  CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
  CPPUNIT_ASSERT(!valid);
}

void LVQDTargetTest::check_orient_invalid( TargetCalculator& Wcalc, 
                                           bool expect_error )
{
  OrientFactor Vcalc(&Wcalc);
  MsqError err;

  MsqMatrix<3,2> V2D;
  bool valid = Vcalc.get_orient_2D( pd2D, 0, Sample(0,0), V2D, err );
  CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
  CPPUNIT_ASSERT(!valid);

  MsqMatrix<3,3> V3D;
  valid = Vcalc.get_orient_3D( pd3D, 0, Sample(0,0), V3D, err );
   CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
 CPPUNIT_ASSERT(!valid);
}

void LVQDTargetTest::check_aspect_invalid( TargetCalculator& Wcalc, 
                                           bool expect_error )

{
  AspectFactor Dcalc(&Wcalc);
  MsqError err;

  MsqMatrix<2,2> D2D;
  bool valid = Dcalc.get_aspect_2D( pd2D, 0, Sample(0,0), D2D, err );
  CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
  CPPUNIT_ASSERT(!valid);

  MsqMatrix<3,3> D3D;
  valid = Dcalc.get_aspect_3D( pd3D, 0, Sample(0,0), D3D, err );
  CPPUNIT_ASSERT_EQUAL(expect_error, err.error());
  CPPUNIT_ASSERT(!valid);
}

void LVQDTargetTest::test_factor_skew_error()
{
  TargetError Wcalc( true );
  check_skew_invalid( Wcalc, true );
}

void LVQDTargetTest::test_factor_size_error()
{
  TargetError Wcalc( true );
  check_size_invalid( Wcalc, true );
}

void LVQDTargetTest::test_factor_orient_error()
{
  TargetError Wcalc( true );
  check_orient_invalid( Wcalc, true );
}

void LVQDTargetTest::test_factor_aspect_error()
{
  TargetError Wcalc( true );
  check_aspect_invalid( Wcalc, true );
}

void LVQDTargetTest::test_LVQD_product_error()
{
  MsqError err;
  bool valid;
  MsqMatrix<3,3> W3D;
  MsqMatrix<3,2> W2D;
  
  TargetError Werr(true);
  ConstantTarget I(1.0);
  SkewFactor Qerr(&Werr), Qvalid(&I);
  SizeFactor Lerr(&Werr), Lvalid(&I);
  OrientFactor Verr(&Werr), Vvalid(&I);
  AspectFactor Derr(&Werr), Dvalid(&I);
  
  LVQDTargetCalculator Linvalid( &Lerr, &Vvalid, &Qvalid, &Dvalid );
  valid = Linvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);
  valid = Linvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);
  
  LVQDTargetCalculator Vinvalid( &Lvalid, &Verr, &Qvalid, &Dvalid );
  valid = Vinvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);
  valid = Vinvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);
  
  LVQDTargetCalculator Qinvalid( &Lvalid, &Vvalid, &Qerr, &Dvalid );
  valid = Qinvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);
  valid = Qinvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);

  LVQDTargetCalculator Dinvalid( &Lvalid, &Vvalid, &Qvalid, &Derr );
  valid = Dinvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);
  valid = Dinvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  CPPUNIT_ASSERT(err.error());
  CPPUNIT_ASSERT(!valid);
}  


void LVQDTargetTest::test_factor_skew_invalid()
{
  TargetError Wcalc( false );
  check_skew_invalid( Wcalc, false );
}

void LVQDTargetTest::test_factor_size_invalid()
{
  TargetError Wcalc( false );
  check_size_invalid( Wcalc, false );
}

void LVQDTargetTest::test_factor_orient_invalid()
{
  TargetError Wcalc( false );
  check_orient_invalid( Wcalc, false );
}

void LVQDTargetTest::test_factor_aspect_invalid()
{
  TargetError Wcalc( false );
  check_aspect_invalid( Wcalc, false );
}

void LVQDTargetTest::test_LVQD_product_invalid()
{
  MsqError err;
  bool valid;
  MsqMatrix<3,3> W3D;
  MsqMatrix<3,2> W2D;
  
  TargetError Werr(false);
  ConstantTarget I(1.0);
  SkewFactor Qerr(&Werr), Qvalid(&I);
  SizeFactor Lerr(&Werr), Lvalid(&I);
  OrientFactor Verr(&Werr), Vvalid(&I);
  AspectFactor Derr(&Werr), Dvalid(&I);
  
  LVQDTargetCalculator Linvalid( &Lerr, &Vvalid, &Qvalid, &Dvalid );
  valid = Linvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
  valid = Linvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
  
  LVQDTargetCalculator Vinvalid( &Lvalid, &Verr, &Qvalid, &Dvalid );
  valid = Vinvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
  valid = Vinvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
  
  LVQDTargetCalculator Qinvalid( &Lvalid, &Vvalid, &Qerr, &Dvalid );
  valid = Qinvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
  valid = Qinvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);

  LVQDTargetCalculator Dinvalid( &Lvalid, &Vvalid, &Qvalid, &Derr );
  valid = Dinvalid.get_2D_target( pd2D, 0, Sample(0,0), W2D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
  valid = Dinvalid.get_3D_target( pd3D, 0, Sample(0,0), W3D, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
}


void LVQDTargetTest::test_factor_skew_zero()
{
  ConstantTarget zero( MsqMatrix<3,3>(0.0), MsqMatrix<3,2>(0.0) );
  check_skew_invalid( zero, false );
}

void LVQDTargetTest::test_factor_orient_zero()
{
  ConstantTarget zero( MsqMatrix<3,3>(0.0), MsqMatrix<3,2>(0.0) );
  check_orient_invalid( zero, false );
}

void LVQDTargetTest::test_factor_aspect_zero()
{
  ConstantTarget zero( MsqMatrix<3,3>(0.0), MsqMatrix<3,2>(0.0) );
  check_aspect_invalid( zero, false );
}


