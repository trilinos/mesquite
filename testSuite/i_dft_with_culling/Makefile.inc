I_DFT_CULLDIR = $(testdir)/i_dft_with_culling

I_DFT_CULLSRC = $(I_DFT_CULLDIR)/main.cpp

I_DFT_CULLEXEC = $(I_DFT_CULLDIR)/main

all_tests: $(I_DFT_CULLEXEC)

$(I_DFT_CULLEXEC): $(I_DFT_CULLSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${I_DFT_CULLSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
