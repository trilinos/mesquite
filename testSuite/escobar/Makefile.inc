ESCOBAR_DIR = $(testdir)/escobar

ESCOBAR_SRC = $(ESCOBAR_DIR)/main.cpp

ESCOBAR_EXEC = $(ESCOBAR_DIR)/main

all_tests: $(ESCOBAR_EXEC)


$(ESCOBAR_EXEC): $(ESCOBAR_SRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) $(CONFIG_LDFLAGS) \
	-o $@ ${ESCOBAR_SRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
