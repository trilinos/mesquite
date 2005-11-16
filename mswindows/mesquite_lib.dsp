# Microsoft Developer Studio Project File - Name="mesquite_lib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=mesquite_lib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mesquite_lib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mesquite_lib.mak" CFG="mesquite_lib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mesquite_lib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "mesquite_lib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "mesquite_lib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "mesquite_lib___Win32_Release"
# PROP BASE Intermediate_Dir "mesquite_lib___Win32_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release_Lib"
# PROP Intermediate_Dir "Release_Lib"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "../include" /I "../src" /I "../src/Control" /I "../src/Control/Wrappers" /I "../src/Mesh" /I "../src/Misc" /I "../src/ObjectiveFunction" /I "../src/QualityAssessor" /I "../src/QualityImprover" /I "../src/QualityImprover/TopologyModifier" /I "../src/QualityImprover/VertexMover" /I "../src/QualityImprover/VertexMover/FeasibleNewton" /I "../src/QualityImprover/VertexMover/ConjugateGradient" /I "../src/QualityImprover/VertexMover/LaplacianSmoothers" /I "../src/QualityImprover/VertexMover/NonSmoothSteepestDescent" /I "../src/QualityImprover/VertexMover/Randomize" /I "../src/QualityImprover/VertexMover/SteepestDescent" /I "../src/QualityMetric" /I "../src/QualityMetric/Shape" /I "../src/QualityMetric/Smoothness" /I "../src/QualityMetric/Untangle" /I "../src/QualityMetric/Volume" /I "../src/TargetCalculator" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "HAVE_CLOCK" /D "HAVE__VSNPRINTF" /D "MESQUITE_STATIC_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Release_Lib\mesquite.lib"

!ELSEIF  "$(CFG)" == "mesquite_lib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "mesquite_lib___Win32_Debug"
# PROP BASE Intermediate_Dir "mesquite_lib___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug_Lib"
# PROP Intermediate_Dir "Debug_Lib"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "../include" /I "../src" /I "../src/Control" /I "../src/Control/Wrappers" /I "../src/Mesh" /I "../src/Misc" /I "../src/ObjectiveFunction" /I "../src/QualityAssessor" /I "../src/QualityImprover" /I "../src/QualityImprover/TopologyModifier" /I "../src/QualityImprover/VertexMover" /I "../src/QualityImprover/VertexMover/FeasibleNewton" /I "../src/QualityImprover/VertexMover/ConjugateGradient" /I "../src/QualityImprover/VertexMover/LaplacianSmoothers" /I "../src/QualityImprover/VertexMover/NonSmoothSteepestDescent" /I "../src/QualityImprover/VertexMover/Randomize" /I "../src/QualityImprover/VertexMover/SteepestDescent" /I "../src/QualityMetric" /I "../src/QualityMetric/Shape" /I "../src/QualityMetric/Smoothness" /I "../src/QualityMetric/Untangle" /I "../src/QualityMetric/Volume" /I "../src/TargetCalculator" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "HAVE_CLOCK" /D "HAVE__VSNPRINTF" /D "MESQUITE_STATIC_LIB" /D "MSQ_TRAP_FPE" /D MSQ_ENABLE_DEBUG="1,2" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Debug_Lib\mesquite.lib"

!ENDIF 

# Begin Target

# Name "mesquite_lib - Win32 Release"
# Name "mesquite_lib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\src\QualityMetric\AddQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\AspectRatioGammaQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFAdd.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFMultiply.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFScalarAdd.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFScalarMultiply.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\ConditionNumberQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\ConjugateGradient\ConjugateGradient.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\CornerTag.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Smoothness\EdgeLengthQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Smoothness\EdgeLengthRangeQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\Exponent.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\FeasibleNewton\FeasibleNewton.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\FileTokenizer.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFT.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\IdealWeightInverseMeanRatio.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\IdealWeightMeanRatio.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Control\InstructionQueue.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\LaplacianSmoothers\LaplacianSmoother.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\LInfTemplate.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Volume\LocalSizeQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\LPtoPTemplate.cpp
# End Source File
# Begin Source File

SOURCE=..\src\TargetCalculator\LVQDTargetCalculator.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\MaxTemplate.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\MeanMidNodeMover.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshImpl.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshImplData.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshImplTags.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MeshTransform.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshWriter.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MesquiteVersion.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqDebug.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqError.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqFPE.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqHessian.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqInterrupt.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MsqMeshEntity.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqTimer.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MsqVertex.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\MultiplyQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\NonSmoothSteepestDescent\NonSmoothSteepestDescent.cpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\ObjectiveFunction.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\PatchData.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\PlanarDomain.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\PowerQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityAssessor\QualityAssessor.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\QualityImprover.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\QualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\Randomize\Randomize.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\RI_DFT.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\ScalarAddQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\ScalarMultiplyQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Control\Wrappers\ShapeImprovementWrapper.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\sI_DFT.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\LaplacianSmoothers\SmartLaplacianSmoother.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\SphericalDomain.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\sRI_DFT.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\SteepestDescent\SteepestDescent.cpp
# End Source File
# Begin Source File

SOURCE=..\src\TargetCalculator\TargetCalculator.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Control\TerminationCriterion.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\TopologyInfo.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\TopologyModifier\TopologyModifier.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Untangle\UntangleBetaQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\Vector3D.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\VertexConditionNumberQualityMetric.cpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\VertexMover.cpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\VtkTypeInfo.cpp
# End Source File
# Begin Source File

SOURCE=..\src\TargetCalculator\WTargetCalculator.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\src\QualityMetric\AddQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\AspectRatioGammaQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFAdd.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFMultiply.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFScalarAdd.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\CompositeOFScalarMultiply.hpp
# End Source File
# Begin Source File

SOURCE=..\src\TargetCalculator\ConcreteTargetCalculators.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\ConditionNumberQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\ConjugateGradient\ConjugateGradient.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\CornerTag.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\DistanceFromTarget.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Smoothness\EdgeLengthQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Smoothness\EdgeLengthRangeQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\Exponent.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\FeasibleNewton\FeasibleNewton.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\FileTokenizer.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFT.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFT_Generalized.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFT_InverseMeanRatio.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFT_NoBarrier.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFT_StrongBarrier.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFT_WeakBarrier.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\I_DFTFamilyFunctions.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\IdealWeightInverseMeanRatio.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\IdealWeightMeanRatio.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Control\InstructionQueue.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Control\Wrappers\LaplacianIQ.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\LaplacianSmoothers\LaplacianSmoother.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\LInfTemplate.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Volume\LocalSizeQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\LPtoPTemplate.hpp
# End Source File
# Begin Source File

SOURCE=..\src\TargetCalculator\LVQDTargetCalculator.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\Matrix3D.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\MaxTemplate.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\MeanMidNodeMover.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\MeanRatioFunctions.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshImpl.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshImplData.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshImplTags.hpp
# End Source File
# Begin Source File

SOURCE=..\include\MeshInterface.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MeshTransform.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MeshWriter.hpp
# End Source File
# Begin Source File

SOURCE=..\include\Mesquite.hpp
# End Source File
# Begin Source File

SOURCE=..\include\mesquite_config.win.h
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqDebug.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqError.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqFPE.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MsqFreeVertexIndexIterator.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqHessian.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqInterrupt.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MsqMeshEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqTag.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\MsqTimer.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\MsqVertex.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\MultiplyQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\NonSmoothSteepestDescent\NonSmoothSteepestDescent.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\NullImprover.hpp
# End Source File
# Begin Source File

SOURCE=..\src\ObjectiveFunction\ObjectiveFunction.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\PatchData.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\PatchDataMem.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\PatchDataUser.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\PlanarDomain.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\PowerQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityAssessor\QualityAssessor.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\QualityImprover.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\QualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\Randomize\Randomize.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\RI_DFT.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\ScalarAddQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\ScalarMultiplyQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Control\Wrappers\ShapeImprovementWrapper.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\ShapeQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\sI_DFT.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\LaplacianSmoothers\SmartLaplacianSmoother.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Smoothness\SmoothnessQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\SphericalDomain.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\DFT\sRI_DFT.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\SteepestDescent\SteepestDescent.hpp
# End Source File
# Begin Source File

SOURCE=..\src\TargetCalculator\TargetCalculator.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\TargetMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Control\TerminationCriterion.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Mesh\TopologyInfo.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\TopologyModifier\TopologyModifier.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Untangle\UntangleBetaQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Untangle\UntangleQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\Vector3D.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Shape\VertexConditionNumberQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityImprover\VertexMover\VertexMover.hpp
# End Source File
# Begin Source File

SOURCE=..\src\QualityMetric\Volume\VolumeQualityMetric.hpp
# End Source File
# Begin Source File

SOURCE=..\src\Misc\VtkTypeInfo.hpp
# End Source File
# Begin Source File

SOURCE=..\src\TargetCalculator\WTargetCalculator.hpp
# End Source File
# End Group
# End Target
# End Project
