# Microsoft Developer Studio Project File - Name="mesquite_dll" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=mesquite_dll - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mesquite_dll.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mesquite_dll.mak" CFG="mesquite_dll - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mesquite_dll - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "mesquite_dll - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "mesquite_dll - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release_DLL"
# PROP BASE Intermediate_Dir "Release_DLL"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release_DLL"
# PROP Intermediate_Dir "Release_DLL"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MESQUITE_DLL_EXPORTS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GR /GX /O2 /I "../include" /I "../src" /I "../src/Control" /I "../src/Control/Wrappers" /I "../src/Mesh" /I "../src/Misc" /I "../src/ObjectiveFunction" /I "../src/QualityAssessor" /I "../src/QualityImprover" /I "../src/QualityImprover/TopologyModifier" /I "../src/QualityImprover/VertexMover" /I "../src/QualityImprover/VertexMover/ConjugateGradient" /I "../src/QualityImprover/VertexMover/FeasibleNewton" /I "../src/QualityImprover/VertexMover/LaplacianSmoothers" /I "../src/QualityImprover/VertexMover/NonSmoothSteepestDescent" /I "../src/QualityImprover/VertexMover/Randomize" /I "../src/QualityImprover/VertexMover/SteepestDescent" /I "../src/QualityMetric" /I "../src/QualityMetric/Shape" /I "../src/QualityMetric/Smoothness" /I "../src/QualityMetric/Untangle" /I "../src/QualityMetric/Volume" /I "../src/TargetCalculator" /D "NDEBUG" /D "USE_STD_INCLUDES" /D "_LIB" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MESQUITE_DLL_EXPORTS" /D "HAVE_CLOCK" /D "HAVE__VSNPRINTF" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386 /out:"Release_DLL/mesquite.dll"

!ELSEIF  "$(CFG)" == "mesquite_dll - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug_DLL"
# PROP BASE Intermediate_Dir "Debug_DLL"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug_DLL"
# PROP Intermediate_Dir "Debug_DLL"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MESQUITE_DLL_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GR /GX /ZI /Od /I "../include" /I "../src" /I "../src/Control" /I "../src/Control/Wrappers" /I "../src/Mesh" /I "../src/Misc" /I "../src/ObjectiveFunction" /I "../src/QualityAssessor" /I "../src/QualityImprover" /I "../src/QualityImprover/TopologyModifier" /I "../src/QualityImprover/VertexMover" /I "../src/QualityImprover/VertexMover/FeasibleNewton" /I "../src/QualityImprover/VertexMover/ConjugateGradient" /I "../src/QualityImprover/VertexMover/LaplacianSmoothers" /I "../src/QualityImprover/VertexMover/NonSmoothSteepestDescent" /I "../src/QualityImprover/VertexMover/Randomize" /I "../src/QualityImprover/VertexMover/SteepestDescent" /I "../src/QualityMetric" /I "../src/QualityMetric/Shape" /I "../src/QualityMetric/Smoothness" /I "../src/QualityMetric/Untangle" /I "../src/QualityMetric/Volume" /I "../src/TargetCalculator" /D "_DEBUG" /D MSQ_ENABLE_DEBUG=1 /D "2" /D "MSQ_TRAP_FPE" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MESQUITE_DLL_EXPORTS" /D "HAVE_CLOCK" /D "HAVE__VSNPRINTF" /D "USE_STD_INCLUDES" /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /out:"Debug_DLL/mesquite.dll" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "mesquite_dll - Win32 Release"
# Name "mesquite_dll - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter ""
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
