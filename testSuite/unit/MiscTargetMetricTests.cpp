#define TARGET_TEST_GROUP "MiscTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "AWSizeNB1.hpp"
#include "AWSizeB1.hpp"
#include "TSizeNB1.hpp"
#include "TSizeB1.hpp"

//                     NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( AWSizeNB1,  true, false,  true, false );
TEST_METRIC_WITH_HESS( AWSizeB1,   true, false,  true,  true );
TEST_METRIC_WITH_HESS( TSizeNB1,  true, false,  true, false );
TEST_METRIC_WITH_GRAD( TSizeB1,   true, false,  true,  true );
