# Change this in Makefile.customize so that it points to your
# installation of cppunit
#CPPUNIT_DIR = ../../../cppunit


# The rest of these items probably won't need to be messed with
#CPPUNIT_INCLUDE = -I${CPPUNIT_DIR}/include
TSTT_INCLUDE = -I../../TSTT-interface
MESQUITE_INCLUDE = ${TSTT_INCLUDE} -I../../include -I../../includeLinks
TEST_INCLUDE = ${CPPUNIT_INCLUDE} ${MESQUITE_INCLUDE}

CPPUNIT_LINK = -L${CPPUNIT_DIR}/src/cppunit/.libs -lcppunit -ldl

MESQUITE_LINK = -L../../lib -lmesquite ${TSTT_LINK}
MESQUITE_AOMD_LINK = -L../../lib -lmesquite ${AOMD_TSTT_LINK}
MESQUITE_MDB_LINK = -L../../lib -lmesquite ${MDB_TSTT_LINK}

TEST_LINK = ${CPPUNIT_LINK} ${MESQUITE_LINK} -lm
TEST_AOMD_LINK = ${CPPUNIT_LINK} ${MESQUITE_AOMD_LINK} -lm
TEST_MDB_LINK = ${CPPUNIT_LINK} ${MESQUITE_MDB_LINK} -lm
PREFIX = @
PREFIX = 

include ../../Makefile.customize

MSQ_BASE_DIR = ../..

TEST_SRCS = msq_test_main.cpp \
    VtkTest.cpp \
    Matrix3DTest.cpp \
    MsqMeshEntityTest.cpp \
    Vector3DTest.cpp \
    InstructionQueueTest.cpp \
    PatchDataTest.cpp \
    QualityMetricTest.cpp \
    ObjectiveFunctionTest.cpp \
    MeshSetTest.cpp \
    MsqVertexTest.cpp \
    MsqFreeVertexIndexIteratorTest.cpp \
    MsqHessianTest.cpp \
    SphericalGeometryTest.cpp \
    PlanarGeometryTest.cpp \
    TerminationCriterionTest.cpp \
    VertexCullingRegressionTest.cpp \
    ExodusTest.cpp \
    MeshInterfaceTest.cpp

TEST_OBJS = ${TEST_SRCS:.cpp=.o}

msq_test : ${TEST_OBJS} ../../lib/libmesquite.a
	@echo "Linking $@"
	${PREFIX}${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_LINK} -o $@

msq_test_aomd : ${TEST_OBJS}
	@echo "Linking $@"
	${PREFIX}${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_AOMD_LINK} -o $@

msq_test_mdb : ${TEST_OBJS}
	@echo "Linking $@"
	${PREFIX}${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_MDB_LINK} -o $@

clean:
	@rm -f *.o msq_test
	@rm -f *.o msq_test_mdb
	@rm -f *.o msq_test_aomd
	@rm -rf SunWS_cache
	@rm -rf cxx_repository

.SUFFIXES : .cpp .o

.cpp.o:
	@echo "Compiling $<"
	${PREFIX}${CXX} -c ${CXXFLAGS} ${TEST_INCLUDE} -o $@ $<
