@UTIL_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la libmsqutil.la
bin_PROGRAMS = vtkrandom vtkxform 
man1_MANS = vtkrandom.1 vtkxform.1 
noinst_LTLIBRARIES = libmsqutil.la

libmsqutil_la_SOURCES = CLArgs.cpp CLArgs.hpp \
                        CLArgFlag.cpp CLArgFlag.hpp \
                        ManPage.cpp ManPage.hpp

vtkrandom_SOURCES = randomize.cpp
vtkxform_SOURCES = transform.cpp

vtkrandom.1: vtkrandom$(EXEEXT)
	./$< -M > $@ || touch $@

vtkxform.1: vtkxform$(EXEEXT)
	./$< -M > $@ || touch $@

MOSTLYCLEANFILES = $(man1_MANS)
