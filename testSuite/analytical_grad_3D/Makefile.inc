ANALYTICALGRAD3DDIR = $(testdir)/analytical_grad_3D

ANALYTICALGRAD3DSRC = $(ANALYTICALGRAD3DDIR)/main.cpp

ANALYTICALGRAD3DEXEC = $(ANALYTICALGRAD3DDIR)/main

all_tests: $(ANALYTICALGRAD3DEXEC)


$(ANALYTICALGRAD3DEXEC): $(ANALYTICALGRAD3DSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${PREFIX} ${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${ANALYTICALGRAD3DSRC} -L$(locallibdir) -lmesquite -L$(AOMD_LIB_DIR) -lAOMD -lm
