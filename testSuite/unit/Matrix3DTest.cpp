// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 12-Nov-02 at 18:05:56
//  LAST-MOD: 19-Dec-02 at 10:08:18 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file Matrix3DTest.cpp

Unit testing of various functions in the Matrix3D class. 

*/
// DESCRIP-END.
//


#include "Vector3D.hpp"
#include "Matrix3D.hpp"

#include <math.h>

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class Matrix3DTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(Matrix3DTest);
  CPPUNIT_TEST (test_equal);
  CPPUNIT_TEST (test_plus);
  CPPUNIT_TEST (test_minus);
  CPPUNIT_TEST (test_transpose);
  CPPUNIT_TEST (test_plus_transpose);
  CPPUNIT_TEST (test_times);
  CPPUNIT_TEST (test_mult_element);
  CPPUNIT_TEST (test_times_vector);
  CPPUNIT_TEST (test_vector_times);
  CPPUNIT_TEST_SUITE_END();

private:
   
  Vector3D e1, e2, e3;
  Matrix3D mIdentity;
  Matrix3D mMat1;
  Matrix3D mMat2;
  Matrix3D mMat2trans;
  Matrix3D mMat1plus2;
  Matrix3D mMat1plus2trans;
  Matrix3D mMat1times2;

public:
  void setUp()
  {
    // sets up the unit vectors
    e1.set(1,0,0);
    e2.set(0,1,0);
    e3.set(0,0,1);
    
    mIdentity =       "1    0    0
                       0    1    0
                       0    0    1";

    mMat1 =           "1    4.2  2
                       5.2  3    4
                       1    7    0.4";
   
    mMat2 =           "2    4    5
                       2    1    3
                       0    7    8";
    
    mMat2trans =      "2    2    0
                       4    1    7
                       5    3    8";
    
    mMat1plus2 =      "3    8.2   7
                       7.2  4     7
                       1    14    8.4";

    mMat1plus2trans = "3    6.2   2
                       9.2  4     11
                       6    10    8.4";

    mMat1times2 =     "10.4 22.2  33.6
                       16.4 51.8  67.0
                       16.0 13.8  29.2";
  }

  void tearDown()
  {
  }
  
public:
  Matrix3DTest()
  {}

  void test_equal()
  {
    CPPUNIT_ASSERT( mMat1==mMat1 );
    CPPUNIT_ASSERT( mMat1!=mMat2 );
  }
   
  void test_plus()
  {
    Matrix3D sum;
    sum = mMat1 + mMat2; 
    CPPUNIT_ASSERT( sum==mMat1plus2 );
  }

  void test_minus()
  {
    Matrix3D res;
    res = mMat1 - mIdentity;
    Matrix3D correct("0    4.2  2
                      5.2  2    4
                      1    7    -0.6");
    CPPUNIT_ASSERT( res==correct );
  }

  void test_transpose()
  {
    Matrix3D trans = transpose(mMat2);
    CPPUNIT_ASSERT( trans==mMat2trans );
  }
   
  void test_plus_transpose()
  {
    Matrix3D plus_trans = mMat1.plus_transpose(mMat2);
    CPPUNIT_ASSERT( plus_trans==mMat1plus2trans );
  }
   
  void test_times()
  {
    Matrix3D mult = mMat1*mMat2;
    CPPUNIT_ASSERT( mult==mMat1times2 );

    mult = mMat1 * mIdentity;
    CPPUNIT_ASSERT( mult==mMat1 );
  }

  void test_mult_element()
  {
    Matrix3D mat = mult_element(mMat1, mIdentity);
    Matrix3D correct("1 0 0
                      0 3 0
                      0 0 0.4");
    CPPUNIT_ASSERT( mat==correct );
  }

  void test_times_vector()
  {
    Vector3D vec = mMat1*e1;
    Vector3D correct(1, 5.2, 1);
    CPPUNIT_ASSERT( vec==correct);
  }

  void test_vector_times()
  {
    Vector3D vec = e1*mMat1;
    Vector3D correct(1, 4.2, 2);
    CPPUNIT_ASSERT( vec==correct );
    
    Vector3D vec2(2.1, 3, 8);
    vec = vec2*mMat1;
    correct.set(25.7, 73.82, 19.4);
    CPPUNIT_ASSERT( vec==correct );
  }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Matrix3DTest, "Matrix3DTest");
 
