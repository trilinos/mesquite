# To switch to AOMD, set this in Makefile.customize:
# TSTT_LINK = -L../../external/AOMD/lib/ -lAOMD -lm

# These are the things you're most likely to need to
# over-ride in Makefile.customize.
CPPUNIT_DIR = ../../../cppunit
MDB_DIR = ../../../mdb
EXODUS_LIB_DIR = ../../../exodus

# The rest of these items probably won't need to be messed with
CPPUNIT_INCLUDE = -I${CPPUNIT_DIR}/include
TSTT_INCLUDE = -I../../TSTT-interface
MESQUITE_INCLUDE = ${TSTT_INCLUDE} -I../../include -I../../includeLinks
INCLUDE = ${CPPUNIT_INCLUDE} ${MESQUITE_INCLUDE}

CPPUNIT_LINK = -L${CPPUNIT_DIR}/src/cppunit/.libs -lcppunit
TSTT_LINK = -L${MDB_DIR} -lMDB -L$(EXODUS_LIB_DIR) -lexoIIv2c -lnetcdf
MESQUITE_LINK = -L../../lib -lmesquite ${TSTT_LINK}
TEST_LINK = ${CPPUNIT_LINK} ${MESQUITE_LINK}

TEST_SRCS = sample_main.cpp \
    MesquiteTestRunner.cpp \
    sample_test.cpp

TEST_OBJS = ${TEST_SRCS:.cpp=.o}

include ../../Makefile.customize

sample_test : ${TEST_OBJS}
	@echo "Linking $@"
	${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_LINK} -o $@

clean:
	@rm -rf *.o sample_test

.SUFFIXES : .cpp .o

.cpp.o:
	@echo "Compiling $<"
	${CXX} -c ${CXXFLAGS} ${TEST_INCLUDE} -o $@ $<
