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
//  LAST-MOD:  5-Dec-02 at 17:09:04 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file ObjectiveFunctionTest.cpp

Unit testing of various functions in the ObjectiveFunction class. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "ObjectiveFunction.hpp"
#include "LPTemplate.hpp"
#include "CompositeOFAdd.hpp"
#include "CompositeOFScalarMultiply.hpp"
#include "GeneralizedConditionNumberQualityMetric.hpp"
#include "MeanRatioQualityMetric.hpp"

#include "PatchDataInstances.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

#include <list>
#include <iterator>

using namespace Mesquite;

class ObjectiveFunctionTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(ObjectiveFunctionTest);
   CPPUNIT_TEST (test_get_quality_metric_list);
   CPPUNIT_TEST_SUITE_END();
   
private:
   
   PatchData m4Quads;

public:
  void setUp()
  {
     MsqError err;
     
     /* our 2D set up: 4 quads, center vertex outcentered by (0,-0.5)
        7____6____5
        |    |    |
        | 2  |  3 |
        8-_  |  _-4       vertex 1 is at (0,0)
        |  -_0_-  |       vertex 5 is at (2,2)
        | 0  |  1 |
        1----2----3
     */
  //   create_four_quads_patch(m4Quads, err); MSQ_CHKERR(err);
  }

  void tearDown()
  {
  }
  
public:
  ObjectiveFunctionTest()
    {}
  
   void test_get_quality_metric_list()
   {
      MsqError err;
      
      // instantiates a couple of QualityMetrics
      ShapeQualityMetric* mean_ratio = MeanRatioQualityMetric::create_new();
      ShapeQualityMetric* condition_nb = GeneralizedConditionNumberQualityMetric::create_new();

      // and creates a composite objective function.
      LPTemplate* LP2_mean_ratio = new LPTemplate(mean_ratio, 2, err); MSQ_CHKERR(err);
      LPTemplate* LP2_condition_nb = new LPTemplate(condition_nb, 2, err); MSQ_CHKERR(err);
      CompositeOFScalarMultiply* LP2_condition_nb_x3 = new CompositeOFScalarMultiply(3,LP2_condition_nb);   
      CompositeOFAdd comp_OF(LP2_mean_ratio, LP2_condition_nb_x3);

      // test the (simple) get_quality_metric function, on OF with 1 QM. 
      QualityMetric* QM_;
      QM_ = LP2_mean_ratio->get_quality_metric();
      CPPUNIT_ASSERT(QM_==mean_ratio);
      // test the (simple) get_quality_metric function, on OF with 2 QMs. 
      QM_ = comp_OF.get_quality_metric();
      CPPUNIT_ASSERT(QM_==NULL);

      // test the get_quality_metric_list function, on OF with 1 QM. 
      std::list<QualityMetric*> QM__;
      QM__ = LP2_mean_ratio->get_quality_metric_list();
      CPPUNIT_ASSERT( QM__.size()==1 );
      CPPUNIT_ASSERT( *(QM__.begin())==mean_ratio );
      // test the get_quality_metric_list function, on OF with 2 QMs. 
      QM__.clear();
      QM__ = comp_OF.get_quality_metric_list();
      CPPUNIT_ASSERT( QM__.size()==2 );
      std::list<QualityMetric*>::const_iterator QM2 = QM__.begin();
      std::list<QualityMetric*>::const_iterator QM1 = QM2++;;
      CPPUNIT_ASSERT( *QM1==mean_ratio   || *QM2==mean_ratio &&
                      *QM1==condition_nb || *QM2==condition_nb  );

      delete mean_ratio;
      delete condition_nb;
      delete LP2_mean_ratio;
      delete LP2_condition_nb;
      delete LP2_condition_nb_x3;
   }
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ObjectiveFunctionTest, "ObjectiveFunctionTest");

