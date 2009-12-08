@TEST_MAKE_INCLUDE@
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@
AM_CPPFLAGS += '-DSRCDIR="$(srcdir)/"'
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt

EXTRA_DIST = expected_linear_output.vtk \
             expected_quadratic_output.vtk \
             hour-quad8.vtk \
	     linear_input.vtk \
	     quadratic_input.vtk

MOSTLYCLEANFILES = smooth_ho.vtk
