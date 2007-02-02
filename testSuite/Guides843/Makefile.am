@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)

MOSTLYCLEANFILES = MeshWithTargetMatrices.vtk \
		   ori_mesh.gpt \
		   ref_mesh.gpt \
		   smo_mesh.vtk \
		   smo_mesh.vtk.gpt
