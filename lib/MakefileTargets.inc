all_libs: $(MSQLIB)

$(MSQLIB): $(MSQLIBSRC)
	@echo "Creating Mesquite Library..."
	@${ARCHIVER} $@ ${MSQLIBSRC}
	@echo "Creation of Mesquite Library complete."


	
