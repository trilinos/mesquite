LAPLWRAPDIR = $(testdir)/laplacian_wrapper_test

LAPLWRAPSRC = $(LAPLWRAPDIR)/main.cpp

LAPLWRAPEXEC = $(LAPLWRAPDIR)/main

all_tests: $(LAPLWRAPEXEC)

$(LAPLWRAPEXEC): $(LAPLWRAPSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${LAPLWRAPSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
