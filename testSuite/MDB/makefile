INCLUDES = -I/disk2/pub/djmelan/mesquite/mdb -I/disk2/pub/djmelan/mesquite/mesquite/includeLinks -I/disk2/pub/djmelan/mesquite/mesquite/include
DEBUG_FLAG = -g
CCFLAGS = -LANG:std ${INCLUDES}
EXODUS_LIB_DIR = /usr/local/eng_sci/cubit/exodus/exodus3.07/lib/irix65_n32
LDFLAGS = -LANG:std -L../../lib -lmesquite -L/disk2/pub/djmelan/mesquite/mdb -lMDB -L$(EXODUS_LIB_DIR) -lexoIIv2c -lnetcdf -lm

mdb_test: main.o
	 ${CXX} ${DEBUG_FLAG} ./main.o ${LDFLAGS} -o $@

.SUFFIXES : .o .cpp

.cpp.o :
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${DEBUG_FLAG} -c ${CCFLAGS} -o $@ $<
