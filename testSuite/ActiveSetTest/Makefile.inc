ACTIVESETTESTDIR = $(testdir)/ActiveSetTest

ACTIVESETTESTSRC = $(ACTIVESETTESTDIR)/main.cpp

ACTIVESETTESTEXEC = $(ACTIVESETTESTDIR)/main

all_tests: $(ACTIVESETTESTEXEC)


$(ACTIVESETTESTEXEC): $(ACTIVESETTESTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${ACTIVESETTESTSRC} -L$(locallibdir) -lmesquite -L$(AOMD_LIB_DIR) -lAOMD -lm
