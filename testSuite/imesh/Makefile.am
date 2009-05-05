@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la $(IMESH_LIBS)
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(top_builddir)/libmesquite.la
TESTS = main

if WITH_IMESH_IMPL

else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	feasible-newton-result.vtk \
			original.vtk \
			smart-laplacian-result.vtk
