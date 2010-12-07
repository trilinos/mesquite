#define TARGET_TEST_GROUP "MiscTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "AWSizeNB1.hpp"
#include "AWSizeB1.hpp"
#include "TSizeNB1.hpp"
#include "TSizeB1.hpp"
#include "TSquared.hpp"

//                     NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( AWSizeNB1,  true, false,  true, false );
TEST_METRIC_WITH_HESS( AWSizeB1,   true, false,  true,  true );
TEST_METRIC_WITH_HESS( TSizeNB1,  true, false,  true, false );
TEST_METRIC_WITH_GRAD( TSizeB1,   true, false,  true,  true );

TSquared test_TSquared;
class TSquared2DTest : public TMetricTest<TSquared,2> {
  public: 
    TSquared2DTest() : TMetricTest<TSquared,2>(false,false,false,false) {}
    CPPUNIT_TEST_SUITE( TSquared2DTest );
    CPPUNIT_TEST( compare_anaytic_and_numeric_grads );
    CPPUNIT_TEST( compare_eval_with_grad_and_eval_with_hess );
    CPPUNIT_TEST( compare_anaytic_and_numeric_hess );
    CPPUNIT_TEST_SUITE_END();
};
class TSquared3DTest : public TMetricTest<TSquared,3> {
  public: 
    TSquared3DTest() : TMetricTest<TSquared,3>(false,false,false,false) {}
    CPPUNIT_TEST_SUITE( TSquared3DTest );
    CPPUNIT_TEST( compare_anaytic_and_numeric_grads );
    CPPUNIT_TEST( compare_eval_with_grad_and_eval_with_hess );
    CPPUNIT_TEST( compare_anaytic_and_numeric_hess );
    CPPUNIT_TEST_SUITE_END();
};
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared3DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, TARGET_TEST_GROUP );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared3DTest, TARGET_TEST_GROUP );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, "TSquaredTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared3DTest, "TSquaredTest" );
