@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = main

MOSTLYCLEANFILES =  ini_mesh.gpt  ref_mesh.gpt  smo_mesh.gpt  smo_mesh.vtk
