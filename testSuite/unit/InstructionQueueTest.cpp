// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD: 13-Nov-02 at 17:31:57 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file InstructionQueueTest.cpp

Unit testing of various functions in the InstructionQueue class. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"
#include "QualityImprover.hpp"
#include "ShapeQualityMetric.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "SteepestDescent.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;

class InstructionQueueTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(InstructionQueueTest);
  CPPUNIT_TEST (test_add_preconditionner);
//  CPPUNIT_TEST ();
  CPPUNIT_TEST_SUITE_END();

private:
   QualityAssessor* mQA;
   QualityImprover* mQI;
   ShapeQualityMetric* mQM;
   ObjectiveFunction* mOF;
   InstructionQueue mQueue;

public:
  void setUp()
  {
     MsqError err;
     // creates a quality assessor and a qualilty improver
     mQM = MeanRatioQualityMetric::create_new();
     mOF = new LPTemplate(mQM, 2, err);
     mQI = new SteepestDescent( mOF );
     mQA = new QualityAssessor(mQM, QualityAssessor::MAXIMUM);
  }

  void tearDown()
  {
     delete mQA;
     delete mQI;
     delete mQM;
     delete mOF;
  }
  
public:
  InstructionQueueTest()
    {}
  
  void test_add_preconditionner()
  {
     MsqError err;
     mQueue.clear();
     mQueue.add_preconditioner(mQI, err);
     CPPUNIT_ASSERT(!err.errorOn);
     mQueue.set_master_quality_improver(mQI,err);
     CPPUNIT_ASSERT(!err.errorOn);
     mQueue.add_preconditioner(mQI, err);
     CPPUNIT_ASSERT_MESSAGE("preconditionner cannot be added after master QI"
                            , err.errorOn);     
  }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(InstructionQueueTest, "InstructionQueueTest");

