/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TestTRel3DMetrics.cpp
 *  \brief Test all implementations of TRel3DMetric
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"

#define TARGET_TEST_GROUP "TestTRel3DMetrics"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "InverseMeanRatio3D.hpp"
#include "TRel3DShape.hpp"
#include "TRel3DShapeBarrier.hpp"
#include "TRel3DShapeBarrierAlt1.hpp"
#include "TRel3DShapeOrient.hpp"
#include "TRel3DShapeOrientAlt1.hpp"
#include "TRel3DShapeOrientBarrier.hpp"
#include "TRel3DShapeOrientBarrierAlt1.hpp"
#include "TRel3DShapeSize.hpp"
#include "TRel3DShapeSizeBarrier.hpp"
#include "TRel3DShapeSizeBarrierAlt1.hpp"
#include "TRel3DShapeSizeBarrierAlt2.hpp"
#include "TRel3DShapeSizeBarrierAlt3.hpp"
#include "TRel3DShapeSizeOrient.hpp"
#include "TRel3DShapeSizeOrientBarrier.hpp"
#include "TRel3DShapeSizeOrientBarrierAlt1.hpp"
#include "TRel3DSize.hpp"
#include "TRel3DSizeBarrier.hpp"
#include "TRel3DUntangleBeta.hpp"
#include "TRel3DUntangleAlt1.hpp"
#include "TRel3DUntangleMu.hpp"
#include "TSquared3D.hpp"

//                     Metric                           !shape !size !orient barrer ideal
TEST_METRIC_WITH_HESS( InverseMeanRatio3D,              false,  true,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DShape,                     false,  true,  true, false );
TEST_METRIC_WITH_HESS( TRel3DShapeBarrier,              false,  true,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeBarrierAlt1,          false,  true,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeOrient,               false,  true, false, false );
TEST_METRIC_WITH_HESS( TRel3DShapeOrientAlt1,           false,  true, false, false );
TEST_METRIC_WITH_HESS( TRel3DShapeOrientBarrier    ,    false,  true, false,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeOrientBarrierAlt1,    false,  true, false,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeSize,                 false, false,  true, false );
TEST_METRIC_WITH_HESS( TRel3DShapeSizeBarrier,          false, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeSizeBarrierAlt1,      false, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeSizeBarrierAlt2,      false, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeSizeBarrierAlt3,      false, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeSizeOrient,           false, false, false, false );
TEST_METRIC_WITH_HESS( TRel3DShapeSizeOrientBarrier,    false, false, false,  true );
TEST_METRIC_WITH_HESS( TRel3DShapeSizeOrientBarrierAlt1,false, false, false,  true );
TEST_METRIC_WITH_HESS( TRel3DSize,                       true, false,  true, false );
TEST_METRIC_WITH_HESS( TRel3DSizeBarrier,                true, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel3DUntangleBeta,               true,  true,  true, false );
TEST_METRIC_WITH_HESS( TRel3DUntangleAlt1,               true,  true,  true, false );
TEST_COMPOSITE_WITH_HESS( TRel3DUntangleMu, TRel3DSize,      true,  false,  true, false );
TEST_COMPOSITE_WITH_HESS( TRel3DUntangleMu, TRel3DShapeSize,false,  false,  true, false );

TSquared3D test_TSquared3D;
class TSquared3DTest : public TMetricTest<TSquared3D> {
  public: 
    TSquared3DTest() : TMetricTest<TSquared3D>(false,false,false,false) {}
    CPPUNIT_TEST_SUITE( TSquared3DTest );
    CPPUNIT_TEST( compare_eval_and_eval_with_grad ); 
    CPPUNIT_TEST( compare_anaytic_and_numeric_grads );
    CPPUNIT_TEST( compare_eval_with_grad_and_eval_with_hess );
    CPPUNIT_TEST( compare_anaytic_and_numeric_hess );
    CPPUNIT_TEST_SUITE_END();
};
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared3DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared3DTest, TARGET_TEST_GROUP );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared3DTest, "TSquared3DTest" );
