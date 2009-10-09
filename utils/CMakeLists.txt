

IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
  INCLUDE(PackageLibraryMacros)
  INCLUDE(PackageAddExecutable)
ENDIF()

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

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

  ADD_LIBRARY( msqutil STATIC ${SOURCES} )
  ADD_EXECUTABLE( vtkrandom randomize.cpp domain.cpp )
  ADD_EXECUTABLE( msqshape  shape-improve.cpp domain.cpp )
  ADD_EXECUTABLE( vtkxform  transform.cpp )
  ADD_EXECUTABLE( msqquality quality.cpp )

ENDIF()

TARGET_LINK_LIBRARIES( vtkrandom  mesquite msqutil )
TARGET_LINK_LIBRARIES( msqshape   mesquite msqutil )
TARGET_LINK_LIBRARIES( vtkxform   mesquite msqutil )
TARGET_LINK_LIBRARIES( msqquality mesquite msqutil )

IF(${PACKAGE_NAME}_ENABLE_SHARED)
  SET( MESQUITE_UTIL_DEFINES "${MESQUITE_DEFINES}" )
ELSE()
  SET( MESQUITE_UTIL_DEFINES "${MESQUITE_DEFINES} -DMESQUITE_STATIC_LIB" )
ENDIF()
SET( UTIL_SOURCES 
     CLArgFlag.cpp 
     CLArgs.cpp 
     domain.cpp 
     ManPage.cpp 
     quality.cpp 
     ramdomize.cpp
     shape-improve.cpp
     transform.cpp )
SET_SOURCE_FILES_PROPERTIES( ${UTIL_SOURCES} COMPILE_FLAGS ${MESQUITE_UTIL_DEFINES} )
