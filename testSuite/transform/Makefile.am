@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)

MOSTLYCLEANFILES = transformed_mesh.vtk
