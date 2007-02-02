@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = main

if WITH_TSTT
  @MESQUITE_TSTTB_MAKE_INCLUDE@
  @MESQUITE_TSTTM_MAKE_INCLUDE@
  LDADD += $(TSTTM_SERVER_LIBS) $(TSTTB_SERVER_LIBS)
else
  TSTT_TEST_SRC =
main$(EXEEXT):
	echo "TSTT support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	feasible-newton-result.vtk \
			original.vtk \
			smart-laplacian-result.vtk
