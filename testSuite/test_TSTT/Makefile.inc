TEST_TSTT_DIR = $(testdir)/test_TSTT

TEST_TSTT_SRC = $(TEST_TSTT_DIR)/main.cpp

TEST_TSTT_EXEC = $(TEST_TSTT_DIR)/main

ifeq ($(MSQ_TSTT_MESH),yes)
all_tests: $(TEST_TSTT_EXEC)
endif

$(TEST_TSTT_EXEC): $(TEST_TSTT_SRC) $(locallibdir)/libmesquite.a \
		   $(TSTT_SHAR_CLIENT_LIB) $(TSTT_MESH_CLIENT_LIB)
	@echo "Linking Mesquite driver application ($@)..."
	@echo "TSTT_BASE_INC=${TSTT_BASE_INC}"
	@echo "TSTT_MESH_INC=${TSTT_MESH_INC}"
	@echo "TSTT_MESH_LNK=${TSTT_MESH_LNK}"
	@echo "TSTT_BASE_LNK=${TSTT_BASE_LNK}"
	@echo "TSTTM_SERVER_LIBS=${TSTTM_SERVER_LIBS}"
	@echo "TSTTB_SERVER_LIBS=${TSTTB_SERVER_LIBS}"
	$(LINKER) $(CXXFLAGS) $(CONFIG_CFLAGS) \
	$(TSTT_BASE_INC) $(TSTT_MESH_INC) -o $@ $< \
	-L$(locallibdir) -lmesquite \
	$(TSTT_MESH_LNK) $(TSTT_BASE_LNK) \
	$(TSTTM_SERVER_LIBS) $(TSTTB_SERVER_LIBS) $(CONFIG_LDFLAGS) -lm
