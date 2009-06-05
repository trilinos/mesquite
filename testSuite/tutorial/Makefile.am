@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = tutorial
tutorial_SOURCES = tutorial.cpp
tutorial_DEPENDENCIES = $(LDADD)

MOSTLYCLEANFILES = original_mesh.vtk  smoothed_mesh.vtk

exampledir = ${docdir}/examples
example_DATA = tutorial.cpp tutorial.make
