all_libs: $(MSQLIB)

$(MSQLIB): $(ALL_MSQ_OBJ)
	@echo "Creating Mesquite Library..."
	$(PREFIX) ${ARCHIVER} $@ $(ALL_MSQ_OBJ) ${EXTRA_MSQ_ARCHIVE_FILES}
	@echo "Creation of Mesquite Library complete."
