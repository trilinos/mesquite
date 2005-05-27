TEST_PYRAMID_DIR = $(testdir)/pyramid

TEST_PYRAMID_SRC = $(TEST_PYRAMID_DIR)/main.cpp

TEST_PYRAMID_EXEC = $(TEST_PYRAMID_DIR)/main

all_tests: $(TEST_PYRAMID_EXEC) 


$(TEST_PYRAMID_EXEC): $(TEST_PYRAMID_SRC) $(locallibdir)/libmesquite.a 
	@echo "Linking Mesquite driver application ($@)..."
	$(LINKER) $(CXXFLAGS) $(LDFLAGS) $(CONFIG_CFLAGS) -o $@ $< \
	-L$(locallibdir) -lmesquite \
	$(CONFIG_LDFLAGS) -lm
