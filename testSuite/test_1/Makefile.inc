TEST1DIR = $(testdir)/test_1

TEST1SRC = $(TEST1DIR)/main.cpp

TEST1EXEC = $(TEST1DIR)/main

all_tests: $(TEST1EXEC)


$(TEST1EXEC): $(TEST1SRC) $(locallibdir)/mesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${TEST1SRC} $(locallibdir)/mesquite.a -Lexternal/AOMD/lib/ -lAOMD -lm
