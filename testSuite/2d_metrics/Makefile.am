@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
AM_CPPFLAGS += '-DSRCDIR="$(srcdir)/"'
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)

EXTRA_DIST = reference.vtk figures.tex
MOSTLYCLEANFILES = *.gpt
