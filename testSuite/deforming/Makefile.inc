all_tests: $(testdir)/deforming/main


$(testdir)/deforming/main: $(testdir)/deforming/main.cpp $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ $< -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
