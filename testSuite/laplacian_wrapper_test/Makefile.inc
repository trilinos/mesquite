LAPLWRAPDIR = $(testdir)/laplacian_wrapper_test

LAPLWRAPSRC = $(LAPLWRAPDIR)/main.cpp

LAPLWRAPEXEC = $(LAPLWRAPDIR)/main

all_tests: $(LAPLWRAPEXEC)

$(LAPLWRAPEXEC): $(LAPLWRAPSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${LAPLWRAPSRC} -L$(locallibdir) -lmesquite -L$(AOMD_LIB_DIR) -lAOMD
