CONVDIR = $(testdir)/convert

CONVV2ESRC = $(CONVDIR)/vtk_to_exo.cpp

CONVV2EEXEC = $(CONVDIR)/vtk_to_exo

CONVE2VSRC = $(CONVDIR)/exo_to_vtk.cpp

CONVE2VEXEC = $(CONVDIR)/exo_to_vtk

all_tests: $(CONVE2VEXEC)
all_tests: $(CONVV2EEXEC)

$(CONVV2EEXEC): $(CONVV2ESRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${CONVV2ESRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
$(CONVE2VEXEC): $(CONVE2VSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${CONVE2VSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
