#define TARGET_TEST_GROUP "UntangleTargetMetricTests"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "TSizeNB1.hpp"
#include "TShapeSize2DNB1.hpp"
#include "TShapeSize3DNB1.hpp"
#include "TMixed.hpp"
#include "TScale.hpp"
#include "TUntangleBeta.hpp"
#include "TUntangle1.hpp"
#include "TUntangleMu.hpp"

class TUntangleShSz : public TMixed
{
public:
  TShapeSize2DNB1 SS2D;
  TShapeSize3DNB1 SS3D;
  TScale SS2DS; // scale 2D value so that it is sensitive to shape deformation
  TUntangleShSz() : TMixed(&SS2DS,&SS3D), SS2DS(10,&SS2D) {}
};

//                               NAME                !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS   ( TUntangleBeta,              true,  true,  true, false );
TEST_METRIC_WITH_HESS   ( TUntangle1,                 true,  true,  true, false );
TEST_COMPOSITE_WITH_HESS( TUntangleMu, TSizeNB1,      true,  false, true, false );
TEST_COMPOSITE_WITH_HESS( TUntangleMu, TUntangleShSz, false, false, true, false );
