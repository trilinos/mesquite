// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Jan-03 at 09:05:56
//  LAST-MOD: 23-Jan-03 at 16:13:56 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqHessianTest.cpp

Unit testing of the MsqHessian class. 

*/
// DESCRIP-END.
//


#include "PatchDataInstances.hpp"
#include "MsqHessian.hpp"

#include <math.h>

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class MsqHessianTest
  : public CppUnit::TestFixture,
    public Mesquite::MsqHessian
{

private:
  CPPUNIT_TEST_SUITE(MsqHessianTest);
  CPPUNIT_TEST (test_initialize);
  CPPUNIT_TEST (test_axpy);
  CPPUNIT_TEST_SUITE_END();

private:
   
  PatchData twoTriangles; 

public:
  void setUp()
  {
    MsqError err;
    create_two_tri_patch(twoTriangles, err); MSQ_CHKERR(err);
  }

  void tearDown()
  {
  }
  
public:
  MsqHessianTest()
  {}

  void test_initialize()
  {
    MsqError err;
    
    MsqHessian::initialize(twoTriangles, err); MSQ_CHKERR(err);

    // Checks values of mRowStart are correct.
    int row_start[] = {0, 4, 7, 8, 9};
    for (int i=0; i<5; ++i) 
      CPPUNIT_ASSERT( mRowStart[i] == row_start[i] );

    // Checks values of mColIndex are correct.
    int col_index[] = {0,1,2,3,1,2,3,2,3};
    for (int i=0; i<9; ++i) 
      CPPUNIT_ASSERT( mColIndex[i] == col_index[i] );

    // Checks values of mAccumulation are correct. 
    int accumulation[] = {0,1,2,4,5,7,0,3,1,8,-6,4};
    for (int i=0; i<12; ++i)
      CPPUNIT_ASSERT( mAccumulation[i] == accumulation[i] );

  }

  void test_accumulate_entries()
  {
    
  }
  
  void test_axpy()
  {
    MsqError err;
    
    MsqHessian::initialize(twoTriangles, err); MSQ_CHKERR(err);

    int hs = MsqHessian::size();

    Vector3D* res = new Vector3D[hs];
    Vector3D* x = new Vector3D[hs];
    Vector3D* y = new Vector3D[hs];

    Matrix3D blocks[6]; // 6 blocks correspond to a triangular element (n+!)n/2 .

    blocks[0] = "4 4 7   4 5 7   7 7 3 ";
    blocks[1] = "4 8 7   3 5 7   1 2 3 ";
    blocks[2] = "4 4 7   6 5 9   1 8 5 ";
    blocks[3] = "4 4 2   4 5 3   2 3 3 ";
    blocks[4] = "2 4 7   3 2 7   1 4 3 ";
    blocks[5] = "8 4 9   4 5 7   9 7 3 ";

    accumulate_entries(twoTriangles, 0, blocks, 6, err); MSQ_CHKERR(err);

    blocks[2] += blocks[5];
    blocks[5] = blocks[3];
    
    accumulate_entries(twoTriangles, 1, blocks, 6, err); MSQ_CHKERR(err);

    Matrix3D entries_6_ans("2 3 1   4 2 4   7 7 3");
    CPPUNIT_ASSERT( mEntries[6] == entries_6_ans );

    x[0].set(4, 5, 6);
    x[1].set(2, 5, 9);
    x[2].set(1, 2, 6);
    x[3].set(1, 5, 9);

    y[0].set(0, 0, 0);
    y[1].set(0, 0, 0);
    y[2].set(0, 0, 0);
    y[3].set(0, 0, 0);

    axpy(res, hs, *this, x, hs, y, hs, err); MSQ_CHKERR(err);

    Vector3D* ans = new Vector3D[hs];
    ans[0].set(636, 635, 453);
    ans[1].set(365, 460, 461);
    ans[2].set(150, 199, 220);
    ans[3].set(166, 204, 174);

    for (int i=0; i<hs; ++i)
      CPPUNIT_ASSERT( res[i] == ans[i] );

    
    y[0].set(3, 2, 6);
    y[1].set(1, 2, 4);
    y[2].set(3, 6, 9);
    y[3].set(2, 4, 4);

    ans[0].set(639, 637, 459);
    ans[1].set(366, 462, 465);
    ans[2].set(153, 205, 229);
    ans[3].set(168, 208, 178);

    axpy(res, hs, *this, x, hs, y, hs, err); MSQ_CHKERR(err);

    for (int i=0; i<hs; ++i)
      CPPUNIT_ASSERT( res[i] == ans[i] );
  }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqHessianTest, "MsqHessianTest");
 
