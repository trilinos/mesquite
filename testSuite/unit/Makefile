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

include ../../Makefile.customize

MSQ_BASE_DIR = ../..

TEST_SRCS = sample_main.cpp \
    MesquiteTestRunner.cpp \
    MsqMeshEntityTest.cpp \
    sample_test.cpp \
    InstructionQueueTest.cpp \
    PatchDataTest.cpp \
    QualityMetricTest.cpp \
    ObjectiveFunctionTest.cpp \
    MeshSetTest.cpp \
    MsqVertexTest.cpp \
    MsqFreeVertexIndexIteratorTest.cpp

TEST_OBJS = ${TEST_SRCS:.cpp=.o}

sample_test : ${TEST_OBJS}
	@echo "Linking $@"
	@${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_LINK} -o $@

clean:
	@rm -rf *.o sample_test
	rm -R SunWS_cache/*

.SUFFIXES : .cpp .o

.cpp.o:
	@echo "Compiling $<"
	@${CXX} -c ${CXXFLAGS} ${TEST_INCLUDE} -o $@ $<
