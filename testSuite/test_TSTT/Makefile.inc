TEST_TSTT_DIR = $(testdir)/test_TSTT

TEST_TSTT_SRC = $(TEST_TSTT_DIR)/main.cpp

TEST_TSTT_EXEC = $(TEST_TSTT_DIR)/main

all_tests: $(TEST_TSTT_EXEC) 


$(TEST_TSTT_EXEC): $(TEST_TSTT_SRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	$(LINKER) $(CXXFLAGS) $(CONFIG_CFLAGS) \
	-o $@ $(TEST_TSTT_SRC) -L$(locallibdir) -lmesquite \
        $(TSTT_LDFLAGS) $(LDFLAGS) $(CONFIG_LDFLAGS) -lm
