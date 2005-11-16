# Microsoft Developer Studio Project File - Name="Unit Tests" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Unit Tests - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "UnitTests.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "UnitTests.mak" CFG="Unit Tests - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Unit Tests - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Unit Tests - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Unit Tests - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GR /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Unit Tests - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "../../../cppunit-1.10.2/include" /I "../include" /I "../src" /I "../src/Control" /I "../src/Control/Wrappers" /I "../src/Mesh" /I "../src/Misc" /I "../src/ObjectiveFunction" /I "../src/QualityAssessor" /I "../src/QualityImprover" /I "../src/QualityImprover/TopologyModifier" /I "../src/QualityImprover/VertexMover" /I "../src/QualityImprover/VertexMover/FeasibleNewton" /I "../src/QualityImprover/VertexMover/ConjugateGradient" /I "../src/QualityImprover/VertexMover/LaplacianSmoothers" /I "../src/QualityImprover/VertexMover/NonSmoothSteepestDescent" /I "../src/QualityImprover/VertexMover/Randomize" /I "../src/QualityImprover/VertexMover/SteepestDescent" /I "../src/QualityMetric" /I "../src/QualityMetric/DFT" /I "../src/QualityMetric/Shape" /I "../src/QualityMetric/Smoothness" /I "../src/QualityMetric/Untangle" /I "../src/QualityMetric/Volume" /I "../src/TargetCalculator" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib mesquite.lib cppunitd.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept /libpath:"Debug" /libpath:"../../../cppunit-1.10.2/lib"
# SUBTRACT LINK32 /nodefaultlib

!ENDIF 

# Begin Target

# Name "Unit Tests - Win32 Release"
# Name "Unit Tests - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\testSuite\unit\ExodusTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\FileTokenizerTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\InstructionQueueTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\Matrix3DTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\MeshInterfaceTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\MeshSetTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\msq_test_main.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\MsqFreeVertexIndexIteratorTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\MsqHessianTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\MsqMeshEntityTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\MsqVertexTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\ObjectiveFunctionTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\PatchDataTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\PlanarGeometryTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\QualityMetricTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\SphericalGeometryTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\TerminationCriterionTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\TopologyInfoTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\Vector3DTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\VertexCullingRegressionTest.cpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\VtkTest.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\testSuite\unit\MesquiteTestRunner.hpp
# End Source File
# Begin Source File

SOURCE=..\testSuite\unit\PatchDataInstances.hpp
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
