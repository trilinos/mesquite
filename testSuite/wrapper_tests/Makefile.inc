WRAPPERTESTSDIR = $(testdir)/wrapper_tests

WRAPPERTESTSSRC = $(WRAPPERTESTSDIR)/main.cpp

WRAPPERTESTSEXEC = $(WRAPPERTESTSDIR)/main

all_tests: $(WRAPPERTESTSEXEC)


$(WRAPPERTESTSEXEC): $(WRAPPERTESTSSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ ${WRAPPERTESTSSRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
