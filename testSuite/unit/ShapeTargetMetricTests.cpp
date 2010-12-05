#define TARGET_TEST_GROUP "ShapeTargetMetricTests"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "TInverseMeanRatio.hpp"
#include "TShape2DNB2.hpp"
#include "TShape3DB2.hpp"
#include "TShapeB1.hpp"
#include "TShapeNB1.hpp"

//                               NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS   ( TInverseMeanRatio, false, true,  true, true  );
TEST_METRIC_WITH_HESS_2D( TShape2DNB2,       false, true,  true, false );
TEST_METRIC_WITH_HESS_3D( TShape3DB2,        false, true,  true, true  );
TEST_METRIC_WITH_HESS   ( TShapeB1,          false, true,  true, true  );
TEST_METRIC_WITH_HESS   ( TShapeNB1,         false, true,  true, false );
