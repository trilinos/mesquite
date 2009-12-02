
IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
  INCLUDE(PackageLibraryMacros)
ENDIF()
  

#
# A) Build w/out support for Fortran name mangling
#

autoconf_header(iBase_FCDefs.h.in iBase_FCDefs.h)
