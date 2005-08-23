TEST_WEDGE_DIR = $(testdir)/wedge

TEST_WEDGE_SRC = $(TEST_WEDGE_DIR)/main.cpp

TEST_WEDGE_EXEC = $(TEST_WEDGE_DIR)/main

all_tests: $(TEST_WEDGE_EXEC) 


$(TEST_WEDGE_EXEC): $(TEST_WEDGE_SRC) $(locallibdir)/libmesquite.a 
	@echo "Linking Mesquite driver application ($@)..."
	$(LINKER) $(CXXFLAGS) $(LDFLAGS) $(CONFIG_CFLAGS) -o $@ $< \
	-L$(locallibdir) -lmesquite \
	$(CONFIG_LDFLAGS) -lm
