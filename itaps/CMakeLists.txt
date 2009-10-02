
IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
  INCLUDE(PackageLibraryMacros)
ENDIF()
  

#
# A) Package-specific configuration options
#

IF(F77_FUNC)
  SET(FC_FUNC 1)
ENDIF()

IF(F77_FUNC_)
  SET(FC_FUNC_ 1)
ENDIF()

autoconf_header(iBase_FCDefs.h.in iBase_FCDefs.h)
