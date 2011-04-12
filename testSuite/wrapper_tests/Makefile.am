@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
AM_CPPFLAGS += "-DSRCDIR=$(srcdir)"
check_PROGRAMS = main deforming
main_SOURCES = main.cpp 
deforming_SOURCES = deforming.cpp
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt sph-10-zsquare.vtk
