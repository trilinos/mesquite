TEST_TSTT_DIR = $(testdir)/test_TSTT

TEST_TSTT_SRC = $(TEST_TSTT_DIR)/main.cpp

TEST_TSTT_EXEC = $(TEST_TSTT_DIR)/main

all_tests: $(TEST_TSTT_EXEC) 


$(TEST_TSTT_EXEC): $(TEST_TSTT_SRC) $(locallibdir)/libmesquite.a \
		   $(TSTT_SHAR_CLIENT_LIB) $(TSTT_MESH_CLIENT_LIB)
	@echo "Linking Mesquite driver application ($@)..."
	$(LINKER) $(CXXFLAGS) $(CONFIG_CFLAGS) -o $@ $< \
	$(TSTT_MESH_CFLAGS) $(TSTT_BASE_INC) $(TSTT_MESH_INC) \
	-L$(locallibdir) -lmesquite $(TSTT_BASE_LNK) $(TSTT_MESH_LNK)\
	$(TSTT_MESH_LDFLAGS) $(CONFIG_LDFLAGS) -lm
