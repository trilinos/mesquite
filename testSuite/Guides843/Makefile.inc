TESTDFTDIR = $(testdir)/Guides843

TESTDFTSRC = $(TESTDFTDIR)/main.cpp

TESTDFTEXEC = $(TESTDFTDIR)/main

all_tests: $(TESTDFTEXEC)


$(TESTDFTEXEC): $(TESTDFTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${TESTDFTSRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
