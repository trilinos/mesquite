all_libs: $(MSQLIB)

$(MSQLIB): $(wildcard $(localobjdir)/*.o )
	@echo "Creating Mesquite Library..."
	@${ARCHIVER} $@ $(wildcard $(localobjdir)/*.o )
	@echo "Creation of Mesquite Library complete."


	
