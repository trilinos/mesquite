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
//  LAST-MOD:  5-Dec-02 at 10:39:44 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqFreeVertexIndexIteratorTest.cpp

Unit testing of various functions in the MsqFreeVertexIndexIterator class. 

 */
// DESCRIP-END.
//


#include "MsqFreeVertexIndexIterator.hpp"
#include "PatchDataInstances.hpp"

#include <math.h>
#include <iostream>

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class MsqFreeVertexIndexIteratorTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqFreeVertexIndexIteratorTest);
  CPPUNIT_TEST (test_hard_fixed_flags);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData pd;

public:
  void setUp()
  {
    MsqError err;

    /*      7____6____5___11
            |    |    |    |
            | 2  |  3 | 5  |
            8-_  |  _-4---10       vertex 1 is at (0,0)
            |  -_0_-  |    |       vertex 11 is at (3,2)
            | 0  |  1 | 4  |
            1----2----3----9
    */
    create_six_quads_patch(pd, err); 
  }

  void tearDown()
  {

  }
  
public:
  MsqFreeVertexIndexIteratorTest()
    {}
  
  void test_hard_fixed_flags()
  {   
     int indices[10];
     int i=0;
     MsqFreeVertexIndexIterator ind(&pd);
     ind.reset();
     while (ind.next()) {
        indices[i] = ind.value();
        ++i;
     } 

     CPPUNIT_ASSERT(indices[0]==0);
     CPPUNIT_ASSERT(indices[1]==4);
     CPPUNIT_ASSERT(i==2); // number of free vertices.

  }


};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqFreeVertexIndexIteratorTest, "MsqFreeVertexIndexIteratorTest");
