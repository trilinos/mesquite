LAPLDIR = $(testdir)/laplacian_test

LAPLSRC = $(LAPLDIR)/main.cpp

LAPLEXEC = $(LAPLDIR)/main

all_tests: $(LAPLEXEC)

$(LAPLEXEC): $(LAPLSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${LAPLSRC} -L$(locallibdir) -lmesquite -L$(AOMD_LIB_DIR) -lAOMD -lm
