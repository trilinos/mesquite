@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/src/libmesquite.la
check_PROGRAMS = par_hex_smooth_laplace
main_SOURCES = par_hex_smooth_laplace.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt

MOSTLYCLEANFILES = *.*.vtk
