default : all

all : mesquite driver

depend : force_rebuild
	@echo "Generating Dependencies for Mesquite Library..."
	@cd src; ${MAKE} depend
	@echo "Creating dependencies in sample driver application..."
	@cd driver; ${MAKE} depend

mesquite : force_rebuild
	@echo "Compiling Mesquite Library..."
	@cd src; ${MAKE}

driver : force_rebuild
	@cd driver; echo "Compiling sample driver application"; ${MAKE}

clean :
	@cd src; ${MAKE} clean
	@cd driver; ${MAKE} clean

force_rebuild :
