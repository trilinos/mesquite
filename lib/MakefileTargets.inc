all_libs: $(MSQLIB)

$(MSQLIB): $(MSQLIBSRC)
	@echo "Creating Mesquite Library..."
	@${ARCHIVER} $@ $(wildcard $(localobjdir)/*.o )
#	@${ARCHIVER} $@ $(MSQLIBSRC)
	@echo "Creation of Mesquite Library complete."


	
