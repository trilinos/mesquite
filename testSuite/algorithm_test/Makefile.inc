ALGTESTDIR = $(testdir)/algorithm_test

ALGTESTSRC = $(ALGTESTDIR)/main.cpp

ALGTESTEXEC = $(ALGTESTDIR)/main

all_tests: $(ALGTESTEXEC)

$(ALGTESTEXEC): $(ALGTESTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${ALGTESTSRC} -L$(locallibdir) -lmesquite -Lexternal/AOMD/lib/ -lAOMD -lm
