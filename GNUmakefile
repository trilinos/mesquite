$$SHELL = /bin/sh

# locallibdir - where the compiled libraries will be linked
locallibdir = ./lib
# locallibdir - where all the include files will be (locally) copied.
localincludedir = ./includeLinks
# localobjdir - where all the object files will be written
localobjdir = ./obj
# file where to write dependencies - will be included in this Makefile
dependenciesfile = make.dependencies
# template_dir - where template instantiations are stored.
#                Only applies to some platforms.
template_dir = 

include Makefile.customize

OUTPUT_OPTION = -o $@ # Use the -o option when generating .o files.
INCLUDE = ${SYSTEM_INCLUDE} -I./include -I$(localincludedir)

# Add to this in each subdirectory.  It's the list of files
# that get added to the mesquite library.
ALL_MSQ_OBJ = 

# Make sure that "default" is the first target
default: depend all

# List all desired module names explicitly below
srcdir = src
MODULENAMES :=  Mesh \
		Control \
		Control/Wrappers \
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
		QualityImprover/VertexMover/FeasibleNewton \
		QualityImprover/VertexMover/Randomize \
		QualityMetric \
		QualityMetric/Shape \
		QualityMetric/Smoothness \
		QualityMetric/Untangle \
		QualityMetric/Volume \
		QualityMetric/DFT \
                TargetCalculator \
		../lib 

testdir = testSuite
TESTNAMES := test_1\
             laplacian_test\
             laplacian_wrapper_test\
             untangle_test\
             analytical_grad_3D\
             algorithm_test\
             ActiveSetTest\
             simple_hybrid_test\
             wrapper_tests\
             tutorial\
             test_DFT \
             test_TSTT \
             escobar \
             Guides843 \
             deforming \
             rezone \
             convert \
             higher_order \
	     transform \
             unit 
###             random_test\
###             convert\
###             shape_wrapper_test

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

include tstt/Makefile.inc

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

settings:
	@echo "TESTS = $(TESTS)"
	@echo "TESTMAKEFILES = $(TESTMAKEFILES)"

depend: 
	@touch $(dependenciesfile)
	@echo "Generating dependencies for all Mesquite source files."
	$(PREFIX) $(MAKEDEPEND)  -f $(dependenciesfile) $(CONFIG_CFLAGS) $(TSTT_INC)\
	$(ALLSRC) 2> /dev/null
	$(PREFIX) cat $(dependenciesfile) | perl -np -e "s/^.*\/(.*\.o:)/obj\/\1/;" > $(dependenciesfile).tmp
	$(PREFIX) mv $(dependenciesfile).tmp $(dependenciesfile)
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
	-rm -f testSuite/unit/*.o testSuite/unit/msq_test
	-rm -f $(TESTNAMES:%=%/main)
	-rm -rf $(template_dir)

veryclean: clean 
	-rm -f $(locallibdir)/*.a $(locallibdir)/*.so
	-rm -f $(dependenciesfile)
	-rm -rf $(localobjdir)/SunWS_cache
	touch $(dependenciesfile)

distrib: all
	@rm -rf mesquite-1.0
	@mkdir mesquite-1.0
	@mkdir mesquite-1.0/lib
	@cp $(MSQLIB) mesquite-1.0/lib
	@mkdir mesquite-1.0/include
	@cp include/*.hpp mesquite-1.0/include
	@cp includeLinks/*.hpp mesquite-1.0/include
	tar cf mesquite-1.0.tar mesquite-1.0
	rm -rf mesquite-1.0

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
