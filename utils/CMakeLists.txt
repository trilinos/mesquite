

IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
  INCLUDE(PackageLibraryMacros)
  INCLUDE(PackageAddExecutable)
ENDIF()

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Mesh)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Misc)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Control)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Control/Wrappers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/MappingFunction)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/MappingFunction/Linear)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/MappingFunction/Lagrange)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityAssessor)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityMetric)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityMetric/Shape)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityMetric/Volume)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityMetric/Untangle)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityMetric/TMP)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityMetric/TMP/2D)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityMetric/TMP/3D)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/TargetCalculator)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/ObjectiveFunction)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityImprover)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityImprover/VertexMover)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityImprover/VertexMover/FeasibleNewton)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityImprover/VertexMover/TrustRegion)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityImprover/VertexMover/ConjugateGradient)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/QualityImprover/VertexMover/Randomize)

APPEND_SET(HEADERS
  CLArgs.hpp
  CLArgFlag.hpp
  ManPage.hpp
  )

APPEND_SET(SOURCES
  CLArgs.cpp
  CLArgFlag.cpp
  ManPage.cpp
  )

IF (TPL_ENABLE_MPI)

ENDIF()

#
# C) Define the targets for package's library(s)
#

IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")

  PACKAGE_ADD_EXECUTABLE(
    vtkrandom
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES randomize.cpp
            domain.hpp
            domain.cpp
    )

  PACKAGE_ADD_EXECUTABLE(
    msqshape
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES shape-improve.cpp
            domain.hpp
            domain.cpp
    )

  PACKAGE_ADD_EXECUTABLE(
    vtkxform
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES transform.cpp
    )

  PACKAGE_ADD_EXECUTABLE(
    msqquality
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES quality.cpp
    )

ELSE()

  ADD_LIBRARY( msqutil ${SOURCES} )
  ADD_EXECUTABLE( vtkrandom randomize.cpp domain.cpp )
  ADD_EXECUTABLE( msqshape  shape-improve.cpp domain.cpp )
  ADD_EXECUTABLE( vtkxform  transform.cpp )
  ADD_EXECUTABLE( msqquality quality.cpp )
  TARGET_LINK_LIBRARIES( vtkrandom  mesquite msqutil )
  TARGET_LINK_LIBRARIES( msqshape   mesquite msqutil )
  TARGET_LINK_LIBRARIES( vtkxform   mesquite msqutil )
  TARGET_LINK_LIBRARIES( msqquality mesquite msqutil )

ENDIF()
