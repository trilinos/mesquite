DEBUG_FLAG = -g
CXX = CC
LINKER = ${CXX}
LDFLAGS = -LANG:std ${DEBUG_FLAG}

CXXFLAGS = -LANG:std ${DEBUG_FLAG}

CPPUNIT_INCLUDE = -I/disk2/pub/djmelan/mesquite/cppunit/include
TSTT_INCLUDE = -I../../TSTT-interface
MESQUITE_INCLUDE = ${TSTT_INCLUDE} -I../../include -I../../includeLinks
TEST_INCLUDE = ${CPPUNIT_INCLUDE} ${MESQUITE_INCLUDE}

CPPUNIT_LINK = -L/disk2/pub/djmelan/mesquite/cppunit/src/cppunit/.libs -lcppunit
EXODUS_LIB_DIR = /usr/local/eng_sci/cubit/exodus/exodus3.07/lib/irix65_n32
TSTT_LINK = -L/disk2/pub/djmelan/mesquite/mdb -lMDB -L$(EXODUS_LIB_DIR) -lexoIIv2c -lnetcdf -lm
MESQUITE_LINK = /disk2/pub/djmelan/mesquite/mesquite/lib/mesquite.a ${TSTT_LINK}
TEST_LINK = ${CPPUNIT_LINK} ${MESQUITE_LINK}

TEST_SRCS = sample_main.cpp \
    MesquiteTestRunner.cpp \
    sample_test.cpp

TEST_OBJS = ${TEST_SRCS:.cpp=.o}

sample_test : ${TEST_OBJS}
	@echo "Linking $@"
	${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_LINK} -o $@

clean:
	@rm -rf *.o sample_test

.SUFFIXES : .cpp .o

.cpp.o:
	@echo "Compiling $<"
	${CXX} -c ${CXXFLAGS} ${TEST_INCLUDE} -o $@ $<
