ACTIVESETTESTDIR = $(testdir)/ActiveSetTest

ACTIVESETTESTSRC = $(ACTIVESETTESTDIR)/main.cpp

ACTIVESETTESTEXEC = $(ACTIVESETTESTDIR)/main

all_tests: $(ACTIVESETTESTEXEC)


$(ACTIVESETTESTEXEC): $(ACTIVESETTESTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ ${ACTIVESETTESTSRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
