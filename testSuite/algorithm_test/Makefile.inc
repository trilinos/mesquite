ALGTESTDIR = $(testdir)/algorithm_test

ALGTESTSRC = $(ALGTESTDIR)/main.cpp

ALGTESTEXEC = $(ALGTESTDIR)/main

all_tests: $(ALGTESTEXEC)

$(ALGTESTEXEC): $(ALGTESTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ ${ALGTESTSRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
