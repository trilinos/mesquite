
all_tests: $(testdir)/rezone/main


$(testdir)/rezone/main: $(testdir)/rezone/main.cpp $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ $< -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
