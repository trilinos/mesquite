#define TARGET_TEST_GROUP "ShapeSizeOrientTargetMetricTests"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "AWShapeSizeOrientNB1.hpp"
#include "TShapeSizeOrientB1.hpp"
#include "TShapeSizeOrientB2.hpp"
#include "TShapeSizeOrientNB1.hpp"

//                            NAME     !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( AWShapeSizeOrientNB1,false,false,false,false );
TEST_METRIC_WITH_HESS( TShapeSizeOrientB1,  false,false,false,true  );
TEST_METRIC_WITH_HESS( TShapeSizeOrientB2,  false,false,false,true  );
TEST_METRIC_WITH_HESS( TShapeSizeOrientNB1, false,false,false,false );
