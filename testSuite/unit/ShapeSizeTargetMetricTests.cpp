#define TARGET_TEST_GROUP "ShapeSizeTargetMetricTests"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "TShapeSize2DB2.hpp"
#include "TShapeSize2DNB1.hpp"
#include "TShapeSize2DNB2.hpp"
#include "TShapeSize3DB2.hpp"
#include "TShapeSize3DB4.hpp"
#include "TShapeSizeB1.hpp"
#include "TShapeSizeB3.hpp"
#include "TShapeSizeNB3.hpp"
#include "TSquared.hpp"

//                            NAME     !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS_2D( TShapeSize2DB2, false,false, true,true  );
TEST_METRIC_WITH_HESS_2D( TShapeSize2DNB1,false,false, true,false );
TEST_METRIC_WITH_HESS_2D( TShapeSize2DNB2,false,false, true,false  );
TEST_METRIC_WITH_HESS_3D( TShapeSize3DB2, false,false, true,true  );
TEST_METRIC_WITH_HESS_3D( TShapeSize3DB4, false,false, true,true  );
TEST_METRIC_WITH_HESS   ( TShapeSizeB1,   false,false, true,true  );
TEST_METRIC_WITH_HESS   ( TShapeSizeB3,   false,false, true,true  );
TEST_METRIC_WITH_HESS   ( TShapeSizeNB3,  false,false, true,false );


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
