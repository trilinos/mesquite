ALGTESTDIR = $(testdir)/algorithm_test

ALGTESTSRC = $(ALGTESTDIR)/main.cpp

ALGTESTEXEC = $(ALGTESTDIR)/main

all_tests: $(ALGTESTEXEC)

$(ALGTESTEXEC): $(ALGTESTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${ALGTESTSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
