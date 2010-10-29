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


/** \file TestTAbs2DMetrics.cpp
 *  \brief Test all implementations of TAbs2DMetric
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"

#define TARGET_TEST_GROUP "TestTAbs2DMetrics"
#include "TMetricTest.hpp"

using namespace Mesquite;

#include "TAbs2DShapeSizeOrient.hpp"

//                     Metric                           !shape !size !orient barrer ideal
TEST_METRIC_WITH_HESS( TAbs2DShapeSizeOrient,           false, false, false, false );
