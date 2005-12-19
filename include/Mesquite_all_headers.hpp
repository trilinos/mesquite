/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 15-Oct-03 at 08:05:56
//  LAST-MOD: 17-Oct-03 at 17:54:14 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*!
  \file   Mesquite_all_headers.hpp
  \brief  This file contains include directives for all mesquite headers.
          Including this header will maximize dependencies and should be
          avoided in the code. Drivers might want to include this header
          to simplify their include directives.

  \author Thomas Leurent
  \date   2003-10-15
*/
// DESCRIP-END.
//

#include "AddQualityMetric.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "BoundedCylinderDomain.hpp"
#include "CompositeOFAdd.hpp"
#include "CompositeOFMultiply.hpp"
#include "CompositeOFScalarAdd.hpp"
#include "CompositeOFScalarMultiply.hpp"
//#include "CompositeQualityMetric.hpp"
#include "ConcreteTargetCalculators.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "ConjugateGradient.hpp"
#include "CylinderDomain.hpp"
#include "DistanceFromTarget.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "EdgeLengthRangeQualityMetric.hpp"
#include "FeasibleNewton.hpp"
#include "GeomTSTT.hpp"
#include "I_DFT.hpp"
#include "I_DFT_Generalized.hpp"
#include "I_DFT_InverseMeanRatio.hpp"
#include "I_DFT_NoBarrier.hpp"
#include "I_DFT_NoBarrierSmoother.hpp"
#include "I_DFT_StrongBarrier.hpp"
#include "I_DFT_WeakBarrier.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "IdealWeightMeanRatio.hpp"
#include "InstructionQueue.hpp"
#include "LInfTemplate.hpp"
#include "LPtoPTemplate.hpp"
#include "LVQDTargetCalculator.hpp"
#include "LaplacianIQ.hpp"
#include "LaplacianSmoother.hpp"
#include "LocalSizeQualityMetric.hpp"
#include "Matrix3D.hpp"
#include "MaxTemplate.hpp"
#include "MeanMidNodeMover.hpp"
#include "MeshImpl.hpp"
#include "MeshTSTT.hpp"
#include "MeshTransform.hpp"
#include "MeshWriter.hpp"
#include "MsqDebug.hpp"
#include "MsqError.hpp"
#include "MsqFPE.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqHessian.hpp"
#include "MsqInterrupt.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqTimer.hpp"
#include "MsqVertex.hpp"
#include "MultiplyQualityMetric.hpp"
#include "NonSmoothSteepestDescent.hpp"
#include "NormSquaredSmoother.hpp"
#include "NullImprover.hpp"
#include "ObjectiveFunction.hpp"
#include "ParameterSet.hpp"
#include "PatchData.hpp"
#include "PatchDataUser.hpp"
#include "PlanarDomain.hpp"
#include "PowerQualityMetric.hpp"
#include "QualityAssessor.hpp"
#include "QualityImprover.hpp"
#include "QualityMetric.hpp"
#include "RI_DFT.hpp"
#include "Randomize.hpp"
#include "ScalarAddQualityMetric.hpp"
#include "ScalarMultiplyQualityMetric.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "ShapeQualityMetric.hpp"
#include "SmartLaplacianSmoother.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "SphericalDomain.hpp"
#include "SteepestDescent.hpp"
#include "TargetCalculator.hpp"
#include "TerminationCriterion.hpp"
#include "TopologyInfo.hpp"
#include "TopologyModifier.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "UntangleQualityMetric.hpp"
#include "Vector3D.hpp"
#include "VertexConditionNumberQualityMetric.hpp"
#include "VertexMover.hpp"
#include "VolumeQualityMetric.hpp"
#include "WTargetCalculator.hpp"
#include "sI_DFT.hpp"
#include "sRI_DFT.hpp"

