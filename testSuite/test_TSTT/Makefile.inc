TESTTSTTDIR = $(testdir)/test_TSTT

TESTTSTTSRC = $(TESTTSTTDIR)/main.cpp

TESTTSTTEXEC = $(TESTTSTTDIR)/main

all_tests: $(TESTTSTTEXEC)


$(TESTTSTTEXEC): $(TESTTSTTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ ${TESTTSTTSRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
