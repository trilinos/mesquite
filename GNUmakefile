$$SHELL = /bin/sh

# locallibdir - where the compiled libraries will be linked
locallibdir = ./lib
# locallibdir - where all the include files will be (locally) copied.
localincludedir = ./includeLinks
# localobjdir - where all the object files will be written
localobjdir = ./obj
# file where to write dependencies - will be included in this Makefile
dependenciesfile = make.dependencies


include Makefile.customize

OUTPUT_OPTION = -o $@ # Use the -o option when generating .o files.
INCLUDE = -I./include -I$(localincludedir) -I./TSTT-interface  

# Make sure that "default" is the first target
default: all depend

# List all desired module names explicitly below
srcdir = src
MODULENAMES :=  Mesh \
		Control \
		Misc \
                ObjectiveFunction \
		QualityAssessor \
		QualityImprover \
		QualityImprover/TopologyModifier \
		QualityImprover/VertexMover \
		QualityImprover/VertexMover/ConjugateGradient \
		QualityImprover/VertexMover/LaplacianSmoothers \
		QualityImprover/VertexMover/NonSmoothSteepestDescent \
		QualityImprover/VertexMover/SteepestDescent \
		QualityImprover/VertexMover/Randomize \
		QualityMetric \
		QualityMetric/Shape \
		QualityMetric/Smoothness \
		QualityMetric/Untangle \
		QualityMetric/Volume \
		../lib 

testdir = testSuite
TESTNAMES := test_1\
             laplacian_test\
             untangle_test\
###             random_test\
###             exo_convert

# ************ inclusion of all the modules 
# ************ MakefileVariables.inc  and
# ************ MakefileTargets.inc

# Add location to the beginning of each (i.e. './')
MODULES := $(wildcard $(patsubst %, $(srcdir)/%, $(MODULENAMES)))

# Generate a list of module makefiles actually present. 
MODULEMAKEFILES := $(wildcard $(patsubst %, %/MakefileVariables.inc,\
                    $(MODULES)))

# the module directories are added to the include path
# not in use - using links setup includeLinks instead 
#INCLUDES += $(patsubst %, -I%, $(MODULES))

# add all .cpp and .cc files to the list of sources. This
# list will be sent to makedepend to automatically generate
# dependancies.
ALLSRC := $(foreach MODULE, $(MODULES),\
	 $(wildcard $(MODULE)/*.cpp $(MODULE)/*.cc)) 

# now include the module makefiles (if there are any)
ifdef MODULEMAKEFILES
include $(MODULEMAKEFILES)
endif

# include any module targets, such as executables, that cannot
# be included in the module makefiles because they require the
# value of variables from other modules (make expands variables
# in the dependancy lists as soon as the target is read)
#MODULETARGETFILES = $(wildcard $(srcdir)/*/MakefileTargets.inc)

MODULETARGETFILES := $(wildcard $(patsubst %, \
	%/MakefileTargets.inc, $(MODULES)))

ifdef MODULETARGETFILES
include $(MODULETARGETFILES)
endif
# ************


# ************ inclusion of all the tests
# ************ Makefile.inc

# Add location to the beginning of each (i.e. './')
TESTS := $(wildcard $(patsubst %, $(testdir)/%, $(TESTNAMES)))

# Generate a list of test makefiles actually present. 
TESTMAKEFILES := $(wildcard $(patsubst %, %/Makefile.inc,\
                    $(TESTS)))

# now include the tests makefiles (if there are any)
ifdef TESTMAKEFILES
include $(TESTMAKEFILES)
endif
# *************


all: all_headers all_objects all_libs 

depend: 
	@touch $(dependenciesfile)
	@echo "Generating dependencies for all Mesquite source files."
	$(MAKEDEPEND)  -f $(dependenciesfile) $(DEPEND_FLAGS) \
	$(ALLSRC) 2> /dev/null
	cat $(dependenciesfile) | perl -np -e "s/^.*\/(.*\.o:)/obj\/\1/;" > $(dependenciesfile).tmp
	mv $(dependenciesfile).tmp $(dependenciesfile)
	@echo " *** Done making depend"


all_headers:
all_objects: all_headers
all_libs: all_objects

#tags:
#	$(ETAGS) `find $(srcdir) -name "*.cc" -o -name "*.hh" \
#		-o -name "*.[chf]"`



clean mostlyclean:
	-rm -f $(foreach MODULE, $(MODULES), $(wildcard $(MODULE)/*.o))
	-rm -f $(localincludedir)/*
	-rm -f $(localobjdir)/*.o

veryclean: clean 
	-rm -f $(locallibdir)/*.a $(locallibdir)/*.so
	-rm -f $(dependenciesfile)
	touch $(dependenciesfile)

#distclean: veryclean
#	-rm -f GNUmakefile config.status config.cache

#GNUmakefile: GNUmakefile.in config.status
#	./config.status

#config.status: configure
#	./config.status --recheck

# This next line is here because of an oddity in the
# SunOS version of 'make'
/opt/SUNWspro/SC5.0/include/CC/rw/ctype :
	@:

include $(dependenciesfile)

# DO NOT DELETE
