HEADERS_DIR = $(testdir)/headers

HEADERS_SRC = $(HEADERS_DIR)/main.cpp

HEADERS_EXEC = $(HEADERS_DIR)/main

all_tests: $(HEADERS_EXEC)


$(HEADERS_EXEC): $(HEADERS_SRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${HEADERS_SRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
