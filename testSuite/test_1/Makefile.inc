TEST1DIR = $(testdir)/test_1

TEST1SRC = $(TEST1DIR)/main.cpp

TEST1EXEC = $(TEST1DIR)/main

all_tests: $(TEST1EXEC)


$(TEST1EXEC): $(TEST1SRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${TEST1SRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
