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

#include "ASMQualityMetric.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "CompositeOFAdd.hpp"
#include "CompositeOFMultiply.hpp"
#include "CompositeOFScalarAdd.hpp"
#include "CompositeOFScalarMultiply.hpp"
#include "CompositeQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "ConjugateGradient.hpp"
#include "CornerJacobianQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "EdgeLengthRangeQualityMetric.hpp"
#include "FeasibleNewton.hpp"
#include "GeneralizedConditionNumberQualityMetric.hpp"
#include "InstructionQueue.hpp"
#include "InverseMeanRatioQualityMetric.hpp"
#include "LaplacianIQ.hpp"
#include "LaplacianSmoother.hpp"
#include "LInfTemplate.hpp"
#include "LocalSizeQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "Matrix3D.hpp"
#include "MaxTemplate.hpp"
#include "MeanRatioFunctions.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "MeshImpl.hpp"
#include "MeshSet.hpp"
#include "MeshTSTT.hpp"
#include "MesquiteError.hpp"
#include "MesquiteInterrupt.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqHessian.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
#include "MsqVertex.hpp"
#include "MultiplyQualityMetric.hpp"
#include "NonSmoothSteepestDescent.hpp"
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
#include "Randomize.hpp"
#include "ScalarAddQualityMetric.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "ShapeQualityMetric.hpp"
#include "SmartLaplacianSmoother.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "SphericalDomain.hpp"
#include "SteepestDescent.hpp"
#include "TerminationCriterion.hpp"
#include "TopologyModifier.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "UntangleQualityMetric.hpp"
#include "Vector3D.hpp"
#include "VertexConditionNumberQualityMetric.hpp"
#include "VertexMover.hpp"
#include "VolumeQualityMetric.hpp"
