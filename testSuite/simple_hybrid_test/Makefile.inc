SIMPLE_HYBRID_DIR = $(testdir)/simple_hybrid_test

SIMPLE_HYBRID_SRC = $(SIMPLE_HYBRID_DIR)/main.cpp

SIMPLE_HYBRID_EXEC = $(SIMPLE_HYBRID_DIR)/main

all_tests: $(SIMPLE_HYBRID_EXEC)


$(SIMPLE_HYBRID_EXEC): $(SIMPLE_HYBRID_SRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${SIMPLE_HYBRID_SRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
