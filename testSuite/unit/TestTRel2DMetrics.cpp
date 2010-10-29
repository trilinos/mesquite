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


/** \file TestTRel2DMetrics.cpp
 *  \brief Test all implementations of TRel2DMetric
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"

#define TARGET_TEST_GROUP "TestTRel2DMetrics"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "InverseMeanRatio2D.hpp"
#include "TRel2DShape.hpp"
#include "TRel2DShapeAlt1.hpp"
#include "TRel2DShapeBarrier.hpp"
#include "TRel2DShapeOrient.hpp"
#include "TRel2DShapeOrientAlt1.hpp"
#include "TRel2DShapeOrientBarrier.hpp"
#include "TRel2DShapeOrientBarrierAlt1.hpp"
#include "TRel2DShapeSize.hpp"
#include "TRel2DShapeSizeAlt1.hpp"
#include "TRel2DShapeSizeAlt2.hpp"
#include "TRel2DShapeSizeBarrier.hpp"
#include "TRel2DShapeSizeBarrierAlt1.hpp"
#include "TRel2DShapeSizeBarrierAlt2.hpp"
#include "TRel2DShapeSizeOrient.hpp"
#include "TRel2DShapeSizeOrientBarrier.hpp"
#include "TRel2DShapeSizeOrientBarrierAlt1.hpp"
#include "TRel2DSize.hpp"
#include "TRel2DSizeBarrier.hpp"
#include "TRel2DUntangleBeta.hpp"
#include "TRel2DUntangleAlt1.hpp"
#include "TRel2DUntangleMu.hpp"
#include "TSquared2D.hpp"

//                     NAME                             !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( InverseMeanRatio2D,              false,  true,  true,  true );
TEST_METRIC_WITH_HESS( TRel2DShape,                     false,  true,  true, false );
TEST_METRIC_WITH_HESS( TRel2DShapeAlt1,                 false,  true,  true, false );
TEST_METRIC_WITH_HESS( TRel2DShapeBarrier,              false,  true,  true,  true );
TEST_METRIC_WITH_HESS( TRel2DShapeOrient,               false,  true, false, false );
TEST_METRIC_WITH_HESS( TRel2DShapeOrientAlt1,           false,  true, false, false );
TEST_METRIC_WITH_HESS( TRel2DShapeOrientBarrier,        false,  true, false,  true );
TEST_METRIC_WITH_HESS( TRel2DShapeOrientBarrierAlt1,    false,  true, false,  true );
TEST_METRIC_WITH_HESS( TRel2DShapeSize,                 false, false,  true, false );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeAlt1,             false, false,  true, false );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeAlt2,             false, false,  true, false );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeBarrier,          false, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeBarrierAlt1,      false, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeBarrierAlt2,      false, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeOrient,           false, false, false, false );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeOrientBarrier,    false, false, false,  true );
TEST_METRIC_WITH_HESS( TRel2DShapeSizeOrientBarrierAlt1,false, false, false,  true );
TEST_METRIC_WITH_HESS( TRel2DSize,                       true, false,  true, false );
TEST_METRIC_WITH_HESS( TRel2DSizeBarrier,                true, false,  true,  true );
TEST_METRIC_WITH_HESS( TRel2DUntangleBeta,               true,  true,  true, false );
TEST_METRIC_WITH_HESS( TRel2DUntangleAlt1,               true,  true,  true, false );
TEST_COMPOSITE_WITH_HESS( TRel2DUntangleMu, TRel2DSize,      true,  false,  true, false );
TEST_COMPOSITE_WITH_HESS( TRel2DUntangleMu, TRel2DShapeSize, true,  false,  true, false );

TSquared2D test_TSquared2D;
class TSquared2DTest : public TMetricTest<TSquared2D> {
  public: 
    TSquared2DTest() : TMetricTest<TSquared2D>(false,false,false,false) {}
    CPPUNIT_TEST_SUITE( TSquared2DTest );
    CPPUNIT_TEST( compare_anaytic_and_numeric_grads );
    CPPUNIT_TEST( compare_eval_with_grad_and_eval_with_hess );
    CPPUNIT_TEST( compare_anaytic_and_numeric_hess );
    CPPUNIT_TEST_SUITE_END();
};
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, "Unit" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, TARGET_TEST_GROUP );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TSquared2DTest, "TSquared2DTest" );
