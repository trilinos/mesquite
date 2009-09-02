
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckIncludeFile)

#
# A) Define the package
#

PACKAGE(Mesquite DISABLE_STRONG_WARNINGS)

#
# B) Set up package-specific options
#

# Read the package version number specified as the second argument
# to the AC_INIT macro in a GNU Autoconf configure.in file.
#
# Define the following variables:
# VERSION_STRING:  The second argument to AC_INIT
# MAJOR_VERSION:   For a version string of the form m.n.p..., m
# MINOR_VERSION:   For a version string of the form m.n.p..., n
# PATCH_VERSION:   For a version string of the form m.n.p..., p

file( READ "configure.in" configure_IN )
string( REGEX REPLACE "^.*AC_INIT *\\([^,]+, *([^, )]+).*$" "\\1" VERSION_STRING "${configure_IN}" )
string( REGEX REPLACE "^([0-9]+)(\\..*)?$" "\\1" MAJOR_VERSION "${VERSION_STRING}" )
string( REGEX REPLACE "^[0-9]+\\.([0-9]+)(\\..*)?$" "\\1" MINOR_VERSION "${VERSION_STRING}" )
if ( VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+" )
  string( REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" PATCH_VERSION "${VERSION_STRING}" )
endif ( VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+" )


SET(MSQ_VERSION ${VERSION_STRING})
SET(MSQ_VERSION_MAJOR ${MAJOR_VERSION})
SET(MSQ_VERSION_MINOR ${MINOR_VERSION})
SET(MSQ_VERSION_PATCH ${PATCH_VERSION})
SET(MSQ_VERSION_STRING "\"Mesquite ${MSQ_VERSION}\"")

#check for specific functions
CHECK_FUNCTION_EXISTS(times     HAVE_TIMES)
CHECK_FUNCTION_EXISTS(clock     HAVE_CLOCK)
CHECK_FUNCTION_EXISTS(cbrt      HAVE_CBRT)
CHECK_FUNCTION_EXISTS(vsnprintf HAVE_VSNPRINTF)
CHECK_FUNCTION_EXISTS(vsprintf  HAVE_VSPRINTF)
CHECK_FUNCTION_EXISTS(finite    HAVE_FINITE)
CHECK_FUNCTION_EXISTS(_Isfinite HAVE__ISFINITE)

#check for header files
CHECK_INCLUDE_FILE(ieeefp.h    HAVE_IEEEFP_H)
CHECK_INCLUDE_FILE(inttypes.h  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE(termios.h   HAVE_TERMIOS_H)
CHECK_INCLUDE_FILE(sys/ioctl.h HAVE_SYS_IOCTL_H)
CHECK_INCLUDE_FILE(dlfcn.h     HAVE_DLFCN_H)
CHECK_INCLUDE_FILE(fcntl.h     HAVE_FCNTL_H)

IF(NOT HAVE_TIMES AND NOT HAVE_CLOCK)
  MESSAGE( FATAL_ERROR
    "\n***"
    "\n*** ERROR: Neither the function times nor the function clock could be found. At least one of these needs to exist or Mesquite cannot build properly."
    "\n***\n"
    )
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_IMESH
  HAVE_Mesquite_IMESH
  "Enable IMESH functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_IGEOM
  HAVE_Mesquite_IGEOM
  "Enable IGEOM functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_IREL
  HAVE_Mesquite_IREL
  "Enable IREL functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_IMESHP
  HAVE_Mesquite_IMESHP
  "Enable IMESHP functionality."
  OFF )

#this needs to be done more like the mesquite version where 
#there is a loop to determine the value. Not sure how to do
#this though. Likely through a try compile.
SET(${PACKAGE_NAME}_MSQ_FUNCTION "" CACHE STRING "Set c++ preprocessor macro for function name.")
SET(MSQ_FUNCTION ${${PACKAGE_NAME}_MSQ_FUNCTION})

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MSQ_USE_OLD_C_HEADERS
  MSQ_USE_OLD_C_HEADERS
  "Enable if old C .h headers must be used for C++."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MSQ_USE_OLD_IO_HEADERS
  MSQ_USE_OLD_IO_HEADERS
  "Enable if old .h headers must be used for C++ I/O."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TRAP_FPE
  MSQ_TRAP_FPE
  "Enable to turn on floating point exception trapping."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_DISABLE_DEBUG_ASSERTS
  NDEBUG
  "Disable asserts in debug mode."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_USE_FUNCTION_TIMERS
  MSQ_USE_FUNCTION_TIMERS
  "Enable function timers."
  OFF )

IF(${PACKAGE_NAME}_ENABLE_DEBUG)
  SET(MSQ_ENABLE_DEBUG)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_ExodusII)
  SET(HAVE_LIBEXOIIV2C ON)
ENDIF()


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(itaps)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(utils)

PACKAGE_ADD_TEST_DIRECTORIES(testSuite)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
