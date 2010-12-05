#define TARGET_TEST_GROUP "UntangleTargetMetricTests"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "TSizeNB1.hpp"
#include "TShapeSizeNB1.hpp"
#include "TUntangleBeta.hpp"
#include "TUntangleAlt1.hpp"
#include "TUntangleMu.hpp"

//                               NAME                !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS   ( TUntangleBeta,              true,  true,  true, false );
TEST_METRIC_WITH_HESS   ( TUntangleAlt1,              true,  true,  true, false );
TEST_COMPOSITE_WITH_HESS( TUntangleMu, TSizeNB1,      true,  false, true, false );
TEST_COMPOSITE_WITH_HESS( TUntangleMu, TShapeSizeNB1, true,  false, true, false );
