/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Target2DTest.cpp
 *  \brief Unit tests for 2D target metrics
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "QualityMetricTester.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "TargetMetric2D.hpp"
#include "JacobianMetric.hpp"
#include "SamplePoints.hpp"
#include "IdealTargetCalculator.hpp"
#include "UnitWeight.hpp"
#include "LinearFunctionSet.hpp"

static const EntityTopology SurfElems[] = { TRIANGLE, QUADRILATERAL };

template <class Metric>
class Target2DTest : public CppUnit::TestFixture
{
private:
  SamplePoints corners;
  LinearFunctionSet mapping;
  QualityMetricTester tester;
  IdealTargetCalculator target;
  UnitWeight weight;
  Metric test_metric;
  JacobianMetric metric;
  bool sizeInvariant, orientInvariant, Barrier;
  double idealVal;
public:
  Target2DTest( bool size_invariant, bool orient_invariant, bool barrier, double ideal_element_val )
    : corners( true, false, false, false ),
      tester( SurfElems, sizeof(SurfElems)/sizeof(SurfElems[0]), &mapping ),
      metric( &corners, &target, &weight, &test_metric, 0 ),
      sizeInvariant(size_invariant), orientInvariant(orient_invariant), Barrier(barrier),
      idealVal(ideal_element_val)
    {}
  
  inline void test_ideal_element_eval() {
    tester.test_evaluate_unit_element( &metric, TRIANGLE, idealVal );
    tester.test_evaluate_unit_element( &metric, QUADRILATERAL, idealVal );
  }
  
  inline void test_ideal_element_gradient() {
    tester.test_ideal_element_zero_gradient( &metric, true );
  }

  inline void test_inverted_element_eval() {
    tester.test_evaluate_inverted_element( &metric, !Barrier );
  }
  
  inline void test_measures_quality() {
    if (sizeInvariant && orientInvariant)
      tester.test_measures_quality( &metric );
  }
  
  inline void test_location_invariant() {
    tester.test_location_invariant( &metric);
    tester.test_grad_location_invariant( &metric );
  }
  
  inline void test_scale() {
    if (sizeInvariant) {
      // these tests is not applicable to the target metrics.
      //tester.test_scale_invariant( &metric );
    }
    else {
      tester.test_measures_size( &metric, true );
    }
  }
  
  inline void test_orient() {
    if (orientInvariant) {
      tester.test_orient_invariant( &metric );
      tester.test_grad_orient_invariant( &metric );
    }
    else {
      tester.test_measures_in_plane_orientation( &metric );
    }
  }
};


#include "Target2DShape.hpp"
#include "Target2DShapeBarrier.hpp"
#include "Target2DShapeOrient.hpp"
#include "Target2DShapeOrientAlt1.hpp"
#include "Target2DShapeOrientAlt2.hpp"
#include "Target2DShapeOrientBarrier.hpp"
#include "Target2DShapeSize.hpp"
#include "Target2DShapeSizeBarrier.hpp"
#include "Target2DShapeSizeBarrierAlt1.hpp"
#include "Target2DShapeSizeBarrierAlt2.hpp"
#include "Target2DShapeSizeOrient.hpp"
#include "Target2DShapeSizeOrientAlt1.hpp"
#include "Target2DShapeSizeOrientBarrier.hpp"
#include "Target2DShapeSizeOrientBarrierAlt2.hpp"

#define REGISTER_TARGET2D_TEST( M, A, B, C, D ) \
class Target2DTest_ ## M : public Target2DTest<M> { public: \
  Target2DTest_ ## M () : Target2DTest<M>( (A), (B), (C), (D) ) {} \
  CPPUNIT_TEST_SUITE( Target2DTest_ ## M ); \
  CPPUNIT_TEST (test_ideal_element_eval); \
  CPPUNIT_TEST (test_ideal_element_gradient); \
  CPPUNIT_TEST (test_inverted_element_eval); \
  CPPUNIT_TEST (test_measures_quality); \
  CPPUNIT_TEST (test_location_invariant); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient); \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Target2DTest_ ## M, "Unit"); 

REGISTER_TARGET2D_TEST( Target2DShape,                      true,  true, false, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeBarrier,               true,  true,  true, 1.0 );
REGISTER_TARGET2D_TEST( Target2DShapeOrient,                true, false, false, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeOrientAlt1,            true, false, false, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeOrientAlt2,            true, false, false, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeOrientBarrier,         true, false,  true, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSize,                 false,  true, false, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSizeBarrier,          false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSizeBarrierAlt1,      false,  true,  true, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSizeBarrierAlt2,      false,  true,  true, 1.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSizeOrient,           false, false, false, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSizeOrientAlt1,       false, false, false, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSizeOrientBarrier,    false, false,  true, 0.0 );
REGISTER_TARGET2D_TEST( Target2DShapeSizeOrientBarrierAlt2,false, false,  true, 0.0 );
