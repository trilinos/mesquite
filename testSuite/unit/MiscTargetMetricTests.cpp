#define TARGET_TEST_GROUP "MiscTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "TSizeNB1.hpp"
#include "TSizeB1.hpp"

//                     NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( TSizeNB1,  true, false,  true, false );
TEST_METRIC_WITH_HESS( TSizeB1,   true, false,  true,  true );
