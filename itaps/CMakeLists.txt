
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

IF(F77_FUNC)
  SET(FC_FUNC 1)
ENDIF()

IF(F77_FUNC_)
  SET(FC_FUNC_ 1)
ENDIF()

PACKAGE_CONFIGURE_FILE(iBase_FCDefs.h)
