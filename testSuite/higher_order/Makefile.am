@TEST_MAKE_INCLUDE@
AM_CPPFLAGS += '-DSRCDIR="$(srcdir)/"'
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
