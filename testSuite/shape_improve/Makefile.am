@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
check_PROGRAMS = main old
main_SOURCES = main.cpp
old_SOURCES = main.cpp
old_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_OLD_WRAPPER
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt
