LAPLDIR = $(testdir)/laplacian_test

LAPLSRC = $(LAPLDIR)/main.cpp

LAPLEXEC = $(LAPLDIR)/main

all_tests: $(LAPLEXEC)

$(LAPLEXEC): $(LAPLSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ ${LAPLSRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
