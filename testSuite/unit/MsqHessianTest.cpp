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
//  LAST-MOD: 14-Jan-03 at 13:32:28 by Thomas Leurent
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
    
    CPPUNIT_ASSERT( mRowStart[0]  == 0 );
    CPPUNIT_ASSERT( mRowStart[1]  == 4 );
    CPPUNIT_ASSERT( mRowStart[2]  == 8 );
    CPPUNIT_ASSERT( mRowStart[3]  == 11 );
    CPPUNIT_ASSERT( mRowStart[4]  == 14 );

    CPPUNIT_ASSERT( mColIndex[0]  == 0 );
    CPPUNIT_ASSERT( mColIndex[1]  == 1 );
    CPPUNIT_ASSERT( mColIndex[2]  == 2 );
    CPPUNIT_ASSERT( mColIndex[3]  == 3 );
    CPPUNIT_ASSERT( mColIndex[4]  == 0 );
    CPPUNIT_ASSERT( mColIndex[5]  == 1 );
    CPPUNIT_ASSERT( mColIndex[6]  == 2 );
    CPPUNIT_ASSERT( mColIndex[7]  == 3 );
    CPPUNIT_ASSERT( mColIndex[8]  == 0 );
    CPPUNIT_ASSERT( mColIndex[9]  == 1 );
    CPPUNIT_ASSERT( mColIndex[10] == 2 );
    CPPUNIT_ASSERT( mColIndex[11] == 0 );
    CPPUNIT_ASSERT( mColIndex[12] == 1 );
    CPPUNIT_ASSERT( mColIndex[13] == 3 );

    CPPUNIT_ASSERT( mColInstr[0]  == 0 );
    CPPUNIT_ASSERT( mColInstr[1]  == 1 );
    CPPUNIT_ASSERT( mColInstr[2]  == 2 );
    CPPUNIT_ASSERT( mColInstr[3]  == 4 );
    CPPUNIT_ASSERT( mColInstr[4]  == 5 );
    CPPUNIT_ASSERT( mColInstr[5]  == 6 );
    CPPUNIT_ASSERT( mColInstr[6]  == 8 );
    CPPUNIT_ASSERT( mColInstr[7]  == 9 );
    CPPUNIT_ASSERT( mColInstr[8]  == 10 );
    CPPUNIT_ASSERT( mColInstr[9]  == 0 );
    CPPUNIT_ASSERT( mColInstr[10] == 3 );
    CPPUNIT_ASSERT( mColInstr[11] == 1 );
    CPPUNIT_ASSERT( mColInstr[12] == 11 );
    CPPUNIT_ASSERT( mColInstr[13] == 13 );
    CPPUNIT_ASSERT( mColInstr[14] == 12 );
    CPPUNIT_ASSERT( mColInstr[15] == 4 );
    CPPUNIT_ASSERT( mColInstr[16] == 7 );
    CPPUNIT_ASSERT( mColInstr[17] == 5 );

  }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqHessianTest, "MsqHessianTest");
 
