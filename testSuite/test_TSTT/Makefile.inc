TESTTSTTDIR = $(testdir)/test_TSTT

TESTTSTTSRC = $(TESTTSTTDIR)/main.cpp

TESTTSTTEXEC = $(TESTTSTTDIR)/main

all_tests: $(TESTTSTTEXEC)


$(TESTTSTTEXEC): $(TESTTSTTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${TESTTSTTSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) $(TSTT_LDFLAGS) -lm
