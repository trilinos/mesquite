UNIT_DIR = $(testdir)/unit

MSQ_TEST = $(UNIT_DIR)/msq_test

all_tests: $(MSQ_TEST)

TEST_SRCS = msq_test_main.cpp \
    FileTokenizerTest.cpp \
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
    MeshInterfaceTest.cpp \
    TopologyInfoTest.cpp 

ifeq ($(MSQ_TSTT_MESH),yes)
    TEST_SRCS += TSTT_Test.cpp
endif

TEST_OBJS = ${TEST_SRCS:.cpp=.o}

TEST_OBJS_P = $(TEST_OBJS:%=$(UNIT_DIR)/%)

$(MSQ_TEST) : $(TEST_OBJS_P) $(locallibdir)/libmesquite.a
	@echo "Linking $@"
	$(PREFIX)$(LINKER) $(LDFLAGS) \
        $(TEST_OBJS_P) \
        -o $@ \
        -L$(locallibdir) -lmesquite \
        $(CPPUNIT_LNK) $(CONFIG_LDFLAGS) $(TSTT_MESH_LDFLAGS)

.SUFFIXES : .cpp .o

.cpp.o:
	@echo "Compiling $<"
	$(PREFIX)$(CXX) -c $(CXXFLAGS) $(CONFIG_CFLAGS) \
        $(MESQUITE_INCLUDE) $(CPPUNIT_INC) $(TSTT_MESH_CFLAGS) \
        -c -o $@ $<
