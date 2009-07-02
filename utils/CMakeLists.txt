
INCLUDE(PackageLibraryMacros)
INCLUDE(PackageAddExecutable)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

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

PACKAGE_ADD_LIBRARY(
  msqutil
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

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

