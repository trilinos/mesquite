# To switch to AOMD, set this in Makefile.customize:
# TSTT_LINK = ${AOMD_TSTT_LINK}

# Change this in Makefile.customize so that it points to your
# installation of cppunit
CPPUNIT_DIR = ../../../cppunit


# The rest of these items probably won't need to be messed with
CPPUNIT_INCLUDE = -I${CPPUNIT_DIR}/include
TSTT_INCLUDE = -I../../TSTT-interface
MESQUITE_INCLUDE = ${TSTT_INCLUDE} -I../../include -I../../includeLinks
INCLUDE = ${CPPUNIT_INCLUDE} ${MESQUITE_INCLUDE}

CPPUNIT_LINK = -L${CPPUNIT_DIR}/src/cppunit/.libs -lcppunit
MESQUITE_LINK = -L../../lib -lmesquite ${TSTT_LINK}
TEST_LINK = ${CPPUNIT_LINK} ${MESQUITE_LINK}
TSTT_LINK = -L../../external/AOMD/lib/ -lAOMD -lm

include ../../Makefile.customize

MSQ_BASE_DIR = ../..

TEST_SRCS = sample_main.cpp \
    MesquiteTestRunner.cpp \
    sample_test.cpp \
    MsqMeshEntityTest.cpp

TEST_OBJS = ${TEST_SRCS:.cpp=.o}

sample_test : ${TEST_OBJS}
	@echo "Linking $@"
	@${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_LINK} ${TSTT_LINK} -o $@

clean:
	@rm -rf *.o sample_test

.SUFFIXES : .cpp .o

.cpp.o:
	@echo "Compiling $<"
	@${CXX} -c ${CXXFLAGS} ${TEST_INCLUDE} -o $@ $<
