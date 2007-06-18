@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = main

if WITH_TEST_IMESH
  @MESQUITE_IMESH_MAKE_INCLUDE@
  LDADD += $(IMESH_CBIND_LIBS_LINK)
else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	feasible-newton-result.vtk \
			original.vtk \
			smart-laplacian-result.vtk
