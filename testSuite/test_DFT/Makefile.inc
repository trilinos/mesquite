TESTDFTDIR = $(testdir)/test_DFT

TESTDFTSRC = $(TESTDFTDIR)/main.cpp

TESTDFTEXEC = $(TESTDFTDIR)/main

all_tests: $(TESTDFTEXEC)


$(TESTDFTEXEC): $(TESTDFTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ ${TESTDFTSRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
