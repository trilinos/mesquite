ALGTESTDIR = $(testdir)/algorithm_test

ALGTESTSRC = $(ALGTESTDIR)/main.cpp

ALGTESTEXEC = $(ALGTESTDIR)/algorithm_test

all_tests: $(ALGTESTEXEC)

$(ALGTESTEXEC): $(ALGTESTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${MYTESTSRC} -L$(locallibdir) -lmesquite -Lexternal/AOMD/lib/ -lAOMD -lm
