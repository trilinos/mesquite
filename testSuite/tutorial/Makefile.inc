TUTORIALDIR = $(testdir)/tutorial

TUTORIALSRC = $(TUTORIALDIR)/tutorial.cpp

TUTORIALEXEC = $(TUTORIALDIR)/tutorial

all_tests: $(TUTORIALEXEC)


$(TUTORIALEXEC): $(TUTORIALSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${TUTORIALSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
