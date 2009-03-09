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

#ifdef MSQ_USE_OLD_IO_HEADERS
# include <iostream.h>
#else
# include <iostream>
#endif

using namespace Mesquite;
const double PI = 3.14159265358979323846;

class LVQDTargetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(LVQDTargetTest);
  CPPUNIT_TEST (test_product_2D);
  CPPUNIT_TEST (test_product_3D);
  CPPUNIT_TEST (test_Q_2D);
  CPPUNIT_TEST (test_Q_3D);
  CPPUNIT_TEST (test_delta_diagonal_2D);
  CPPUNIT_TEST (test_delta_diagonal_3D);
  CPPUNIT_TEST (test_scale_2D);
  CPPUNIT_TEST (test_scale_3D);
  CPPUNIT_TEST (test_orient_2D);
  CPPUNIT_TEST (test_orient_3D);
  CPPUNIT_TEST (test_shape_2D);
  CPPUNIT_TEST (test_shape_3D);
  CPPUNIT_TEST (test_aspect_ratio_2D);
  CPPUNIT_TEST (test_aspect_ratio_3D);
  CPPUNIT_TEST_SUITE_END();

public:

  void test_product_2D();
  void test_product_3D();
  void test_Q_2D();
  void test_Q_3D();
  void test_delta_diagonal_2D();
  void test_delta_diagonal_3D();
  void test_scale_2D();
  void test_scale_3D();
  void test_orient_2D();
  void test_orient_3D();
  void test_shape_2D();
  void test_shape_3D();
  void test_aspect_ratio_2D();
  void test_aspect_ratio_3D();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LVQDTargetTest, "LVQDTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LVQDTargetTest, "Unit");

class TargetCalculator3D : public TargetCalculator
{
  private:
    MsqMatrix<3,3> target;
    bool surf3D;
  
  public:
    TargetCalculator3D( const MsqMatrix<3,3>& value, bool surf_3d = false )
      : target(value), surf3D(surf_3d) {}
    
    virtual bool get_3D_target( PatchData&, size_t, Sample, MsqMatrix<3,3>& result, MsqError& )
      { result = target; return true; }
    
    virtual bool get_2D_target( PatchData&, size_t, Sample, MsqMatrix<3,2>& result, MsqError& )
      { CPPUNIT_ASSERT(false); return false; }
      
    virtual bool surface_targets_are_3D() const
      { return surf3D; }
};

class TargetCalculator2D : public TargetCalculator
{
  private:
    MsqMatrix<3,2> target;
  
  public:
    TargetCalculator2D( const MsqMatrix<3,2>& value ) : target(value) {}
    
    virtual bool get_3D_target( PatchData&, size_t, Sample, MsqMatrix<3,3>& result, MsqError& )
      { CPPUNIT_ASSERT(false); return false; }
    
    virtual bool get_2D_target( PatchData&, size_t, Sample, MsqMatrix<3,2>& result, MsqError& )
      { result = target; return true; }
      
    virtual bool surface_targets_are_3D() const
      { return false; }
};


void LVQDTargetTest::test_product_2D()
{
  MsqPrintError err( msq_stdio::cout );
  const double m2[] = { 1, 2, 
                       -1, 0,
                        1, 1 }; 
  MsqMatrix<3,2> t1(1.0), t2(m2), W;
  PatchData pd;
  bool rval;
  
  TargetCalculator2D tc1(t1);
  LVQDTargetCalculator lvqd1( &tc1, &tc1, &tc1, &tc1 );
  rval = lvqd1.get_2D_target( pd, 0, Sample(0,0), W, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  ASSERT_MATRICES_EQUAL( t1, W, 1e-6 );
  
  TargetCalculator2D tc2(t2);
  LVQDTargetCalculator lvqd2( &tc2, &tc2, &tc2, &tc2 );
  rval = lvqd2.get_2D_target( pd, 0, Sample(0,0), W, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  ASSERT_MATRICES_EQUAL( t2, W, 1e-6 );
}


void LVQDTargetTest::test_product_3D()
{
  MsqPrintError err( msq_stdio::cout );
  const double m2[] = { 1, 2, -1,
                       -1, 0,  1,
                        1, 1,  3 }; 
  MsqMatrix<3,3> t1(1.0), t2(m2), W;
  PatchData pd;
  bool rval;
  
  TargetCalculator3D tc1(t1);
  LVQDTargetCalculator lvqd1( &tc1, &tc1, &tc1, &tc1 );
  rval = lvqd1.get_3D_target( pd, 0, Sample(0,0), W, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  ASSERT_MATRICES_EQUAL( t1, W, 1e-6 );
  
  TargetCalculator3D tc2(t2);
  LVQDTargetCalculator lvqd2( &tc2, &tc2, &tc2, &tc2 );
  rval = lvqd2.get_3D_target( pd, 0, Sample(0,0), W, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  ASSERT_MATRICES_EQUAL( t2, W, 1e-6 );
}

void LVQDTargetTest::test_Q_2D()
{
  const double m2[] = { 1, 2, 
                       -1, 0,
                        1, 1 }; 
  MsqMatrix<3,2> t1(1.0), t2(m2);
  MsqMatrix<2,2> Q;
  
  Q = LVQDTargetCalculator::calc_Q_2D( t1 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(1,0), DBL_EPSILON );
  
  Q = LVQDTargetCalculator::calc_Q_2D( t2 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(1,0), DBL_EPSILON );
}

void LVQDTargetTest::test_Q_3D()
{
  const double m2[] = { 1, 2, -1,
                       -1, 0,  1,
                        1, 1,  3 }; 
  MsqMatrix<3,3> t1(1.0), t2(m2), Q;
  
  Q = LVQDTargetCalculator::calc_Q_3D( t1 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(2,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(2,1), DBL_EPSILON );
  
  Q = LVQDTargetCalculator::calc_Q_3D( t2 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(2,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(2,1), DBL_EPSILON );
}

void LVQDTargetTest::test_delta_diagonal_2D()
{
  const double m2[] = { 1, 2, 
                       -1, 0,
                        1, 1 }; 
  MsqMatrix<3,2> t1(1.0), t2(m2);
  MsqMatrix<2,2> D;
  
  D = LVQDTargetCalculator::calc_delta_2D( t1 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT( sqr_Frobenius(D) > DBL_EPSILON );
  
  D = LVQDTargetCalculator::calc_delta_2D( t2 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT( sqr_Frobenius(D) > DBL_EPSILON );
}

void LVQDTargetTest::test_delta_diagonal_3D()
{
  const double m2[] = { 1, 2, -1,
                       -1, 0,  1,
                        1, 1,  3 }; 
  MsqMatrix<3,3> t1(1.0), t2(m2), D;
  
  D = LVQDTargetCalculator::calc_delta_3D( t1 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(0,2), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(1,2), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(2,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(2,1), DBL_EPSILON );
  CPPUNIT_ASSERT( sqr_Frobenius(D) > DBL_EPSILON );
  
  D = LVQDTargetCalculator::calc_delta_3D( t2 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(0,2), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(1,2), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(2,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, D(2,1), DBL_EPSILON );
  CPPUNIT_ASSERT( sqr_Frobenius(D) > DBL_EPSILON );
}

void LVQDTargetTest::test_scale_2D()
{
  const double m2[] = { 1, 2, 
                       -1, 0,
                        1, 1 }; 
  MsqMatrix<3,2> t1(1.0), t2(m2), V1, V2, S;
  MsqMatrix<2,2> Q1, Q2, D1, D2;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_2D( t1 );
  V1 = LVQDTargetCalculator::calc_V_2D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_2D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_2D( t1 );
  S = 2.0 * t1;
  l2 = LVQDTargetCalculator::calc_lambda_2D( S );
  V2 = LVQDTargetCalculator::calc_V_2D( S );
  Q2 = LVQDTargetCalculator::calc_Q_2D( S );
  D2 = LVQDTargetCalculator::calc_delta_2D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0*l1, l2, 1e-6 );
  ASSERT_MATRICES_EQUAL( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
  
  l1 = LVQDTargetCalculator::calc_lambda_2D( t2 );
  V1 = LVQDTargetCalculator::calc_V_2D( t2 );
  Q1 = LVQDTargetCalculator::calc_Q_2D( t2 );
  D1 = LVQDTargetCalculator::calc_delta_2D( t2 );
  S = 0.1 * t2;
  l2 = LVQDTargetCalculator::calc_lambda_2D( S );
  V2 = LVQDTargetCalculator::calc_V_2D( S );
  Q2 = LVQDTargetCalculator::calc_Q_2D( S );
  D2 = LVQDTargetCalculator::calc_delta_2D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1*l1, l2, 1e-6 );
  ASSERT_MATRICES_EQUAL( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
}

void LVQDTargetTest::test_scale_3D()
{
  const double m2[] = { 1, 2, -1,
                       -1, 0,  1,
                        1, 1,  3 }; 
  MsqMatrix<3,3> t1(1.0), t2(m2), V1, V2, Q1, Q2, D1, D2, S;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_3D( t1 );
  V1 = LVQDTargetCalculator::calc_V_3D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_3D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_3D( t1 );
  S = 2.0 * t1;
  l2 = LVQDTargetCalculator::calc_lambda_3D( S );
  V2 = LVQDTargetCalculator::calc_V_3D( S );
  Q2 = LVQDTargetCalculator::calc_Q_3D( S );
  D2 = LVQDTargetCalculator::calc_delta_3D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0*l1, l2, 1e-6 );
  ASSERT_MATRICES_EQUAL( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
  
  l1 = LVQDTargetCalculator::calc_lambda_3D( t2 );
  V1 = LVQDTargetCalculator::calc_V_3D( t2 );
  Q1 = LVQDTargetCalculator::calc_Q_3D( t2 );
  D1 = LVQDTargetCalculator::calc_delta_3D( t2 );
  S = 0.1 * t2;
  l2 = LVQDTargetCalculator::calc_lambda_3D( S );
  V2 = LVQDTargetCalculator::calc_V_3D( S );
  Q2 = LVQDTargetCalculator::calc_Q_3D( S );
  D2 = LVQDTargetCalculator::calc_delta_3D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1*l1, l2, 1e-6 );
  ASSERT_MATRICES_EQUAL( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
}

static MsqMatrix<3,3> rotation( double a, double x, double y, double z )
{
  MsqMatrix<3,1> axis;
  axis(0,0) = x; axis(1,0) = y; axis(2,0) = z;
  axis /= length(axis);
  
  const double c = cos(a);
  const double s = sin(a);
  MsqMatrix<3,3> R(c);
  R(0,1) = -axis(2,0) * s;
  R(0,2) =  axis(1,0) * s;
  R(1,2) = -axis(0,0) * s;
  R(1,0) =  axis(2,0) * s;
  R(2,0) = -axis(1,0) * s;
  R(2,1) =  axis(0,0) * s;
  return R + (1.0 - c)*outer( axis, axis );
}

void LVQDTargetTest::test_orient_2D()
{
  const double m2[] = { 1, 2, 
                       -1, 0,
                        1, 1 }; 
  MsqMatrix<3,2> t1(1.0), t2(m2), V1, V2, S;
  MsqMatrix<2,2> Q1, Q2, D1, D2;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_2D( t1 );
  V1 = LVQDTargetCalculator::calc_V_2D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_2D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_2D( t1 );
  S = rotation(PI/4, 0, 0, 1) * t1;
  l2 = LVQDTargetCalculator::calc_lambda_2D( S );
  V2 = LVQDTargetCalculator::calc_V_2D( S );
  Q2 = LVQDTargetCalculator::calc_Q_2D( S );
  D2 = LVQDTargetCalculator::calc_delta_2D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( l1, l2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
  
  l1 = LVQDTargetCalculator::calc_lambda_2D( t2 );
  V1 = LVQDTargetCalculator::calc_V_2D( t2 );
  Q1 = LVQDTargetCalculator::calc_Q_2D( t2 );
  D1 = LVQDTargetCalculator::calc_delta_2D( t2 );
  S = rotation(-PI/3, 0, 0, 1) * t2;
  l2 = LVQDTargetCalculator::calc_lambda_2D( S );
  V2 = LVQDTargetCalculator::calc_V_2D( S );
  Q2 = LVQDTargetCalculator::calc_Q_2D( S );
  D2 = LVQDTargetCalculator::calc_delta_2D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( l1, l2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
}



void LVQDTargetTest::test_orient_3D()
{
  const double m2[] = { 1, 2, -1,
                       -1, 0,  1,
                        1, 1,  3 }; 
  MsqMatrix<3,3> t1(1.0), t2(m2), V1, V2, Q1, Q2, D1, D2, S;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_3D( t1 );
  V1 = LVQDTargetCalculator::calc_V_3D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_3D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_3D( t1 );
  S = rotation( PI/4, 1, 1, 1 ) * t1;
  l2 = LVQDTargetCalculator::calc_lambda_3D( S );
  V2 = LVQDTargetCalculator::calc_V_3D( S );
  Q2 = LVQDTargetCalculator::calc_Q_3D( S );
  D2 = LVQDTargetCalculator::calc_delta_3D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( l1, l2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
  
  l1 = LVQDTargetCalculator::calc_lambda_3D( t2 );
  V1 = LVQDTargetCalculator::calc_V_3D( t2 );
  Q1 = LVQDTargetCalculator::calc_Q_3D( t2 );
  D1 = LVQDTargetCalculator::calc_delta_3D( t2 );
  S = rotation( PI/3, 1, 0, 0 ) * t2;
  l2 = LVQDTargetCalculator::calc_lambda_3D( S );
  V2 = LVQDTargetCalculator::calc_V_3D( S );
  Q2 = LVQDTargetCalculator::calc_Q_3D( S );
  D2 = LVQDTargetCalculator::calc_delta_3D( S );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( l1, l2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
}

void LVQDTargetTest::test_shape_2D()
{
  const double m1[] = { 2, 2, 
                       -2, 0,
                        1, 1 }; 
  const double m2[] = { 0, 2,
                       -3, 0,
                        0, 1 };
  MsqMatrix<3,2> t1(m1), t2(m2), V1, V2;
  MsqMatrix<2,2> Q1, Q2, D1, D2;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_2D( t1 );
  V1 = LVQDTargetCalculator::calc_V_2D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_2D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_2D( t1 );
  l2 = LVQDTargetCalculator::calc_lambda_2D( t2 );
  V2 = LVQDTargetCalculator::calc_V_2D( t2 );
  Q2 = LVQDTargetCalculator::calc_Q_2D( t2 );
  D2 = LVQDTargetCalculator::calc_delta_2D( t2 );
  CPPUNIT_ASSERT( l1 < l2 );
  ASSERT_MATRICES_DIFFERENT( V1, V2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
}

void LVQDTargetTest::test_shape_3D()
{
  const double m1[] = { 2, 2, -1.0/sqrt(5.0),
                       -2, 0, 0,
                        1, 1, 2.0/sqrt(5.0) }; 
  const double m2[] = { 0, 2, -1.0/sqrt(5.0),
                       -3, 0, 0,
                        0, 1, 2.0/sqrt(5.0) };
  MsqMatrix<3,3> t1(m1), t2(m2), V1, V2, Q1, Q2, D1, D2;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_3D( t1 );
  V1 = LVQDTargetCalculator::calc_V_3D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_3D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_3D( t1 );
  l2 = LVQDTargetCalculator::calc_lambda_3D( t2 );
  V2 = LVQDTargetCalculator::calc_V_3D( t2 );
  Q2 = LVQDTargetCalculator::calc_Q_3D( t2 );
  D2 = LVQDTargetCalculator::calc_delta_3D( t2 );
  CPPUNIT_ASSERT( l1 < l2 );
  ASSERT_MATRICES_DIFFERENT( V1, V2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_EQUAL( D1, D2, 1e-6 );
}

void LVQDTargetTest::test_aspect_ratio_2D()
{
  const double m1[] = { 2, 2, 
                       -2, 0,
                        1, 1 }; 
  const double m2[] = { 1, 2,
                       -1, 0,
                      0.5, 1 };
  MsqMatrix<3,2> t1(m1), t2(m2), V1, V2;
  MsqMatrix<2,2> Q1, Q2, D1, D2;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_2D( t1 );
  V1 = LVQDTargetCalculator::calc_V_2D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_2D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_2D( t1 );
  l2 = LVQDTargetCalculator::calc_lambda_2D( t2 );
  V2 = LVQDTargetCalculator::calc_V_2D( t2 );
  Q2 = LVQDTargetCalculator::calc_Q_2D( t2 );
  D2 = LVQDTargetCalculator::calc_delta_2D( t2 );
  CPPUNIT_ASSERT( l1 > l2 );
  ASSERT_MATRICES_EQUAL( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( D1, D2, 1e-6 );
}

void LVQDTargetTest::test_aspect_ratio_3D()
{
  const double m1[] = { 1, 2, -1,
                       -1, 0,  1,
                        1, 1,  3 }; 
  const double m2[] = { 1, 1, -1,
                       -1, 0,  1,
                        1,0.5, 3 }; 
  MsqMatrix<3,3> t1(m1), t2(m2), V1, V2, Q1, Q2, D1, D2;
  double l1, l2;
  
  l1 = LVQDTargetCalculator::calc_lambda_3D( t1 );
  V1 = LVQDTargetCalculator::calc_V_3D( t1 );
  Q1 = LVQDTargetCalculator::calc_Q_3D( t1 );
  D1 = LVQDTargetCalculator::calc_delta_3D( t1 );
  l2 = LVQDTargetCalculator::calc_lambda_3D( t2 );
  V2 = LVQDTargetCalculator::calc_V_3D( t2 );
  Q2 = LVQDTargetCalculator::calc_Q_3D( t2 );
  D2 = LVQDTargetCalculator::calc_delta_3D( t2 );
  CPPUNIT_ASSERT( l1 > l2 );
  ASSERT_MATRICES_EQUAL( V1, V2, 1e-6 );
  ASSERT_MATRICES_EQUAL( Q1, Q2, 1e-6 );
  ASSERT_MATRICES_DIFFERENT( D1, D2, 1e-6 );
}
  

