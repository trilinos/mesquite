ACTIVESETTESTDIR = $(testdir)/ActiveSetTest

ACTIVESETTESTSRC = $(ACTIVESETTESTDIR)/main.cpp

ACTIVESETTESTEXEC = $(ACTIVESETTESTDIR)/main

all_tests: $(ACTIVESETTESTEXEC)


$(ACTIVESETTESTEXEC): $(ACTIVESETTESTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${ACTIVESETTESTSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
