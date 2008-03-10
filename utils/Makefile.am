@UTIL_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la libmsqutil.la
bin_PROGRAMS = vtkrandom vtkxform #vtkskin vtk2exo exo2vtk
man1_MANS = vtkrandom.1 vtkxform.1 #vtkskin.man vtk2exo.man exo2vtk.man
noinst_LTLIBRARIES = libmsqutil.la

libmsqutil_la_SOURCES = CLArgs.cpp CLArgs.hpp \
                        CLArgFlag.cpp CLArgFlag.hpp \
                        ManPage.cpp ManPage.hpp
vtkrandom_SOURCES = randomize.cpp
vtkxform_SOURCES = transform.cpp
#vtkskin_SOURCES = skin.cpp
#vtk2exo_SOURCES = vtk2exo.cpp
#exo2vtk_SOURCES = exo2vtk.cpp

vtkrandom.1: vtkrandom
	$< -M > $@

vtkxform.1: vtkxform
	$< -M > $@

vtkskin.1: vtkskin
	$< -M > $@

vtk2exo.1: vtk2exo
	$< -M > $@

exo2vtk.1: exo2vtk
	$< -M > $@
