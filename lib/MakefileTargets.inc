all_libs: $(MSQLIB)

$(MSQLIB): $(ALL_MSQ_OBJ)
	@echo "Creating Mesquite Library..."
	@${ARCHIVER} $@ $(ALL_MSQ_OBJ)
	@echo "Creation of Mesquite Library complete."
