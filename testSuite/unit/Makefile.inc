UNIT_DIR = $(testdir)/unit

MSQ_TEST = $(UNIT_DIR)/msq_test

all_tests: $(MSQ_TEST)

TEST_SRCS = msq_test_main.cpp \
    ExodusTest.cpp \
    FileTokenizerTest.cpp \
    InstructionQueueTest.cpp \
    Matrix3DTest.cpp \
    MeshInterfaceTest.cpp \
    MeshSetTest.cpp \
    MsqFreeVertexIndexIteratorTest.cpp \
    MsqHessianTest.cpp \
    MsqMeshEntityTest.cpp \
    MsqVertexTest.cpp \
    ObjectiveFunctionTest.cpp \
    PatchDataTest.cpp \
    PlanarGeometryTest.cpp \
    QualityMetricTest.cpp \
    SphericalGeometryTest.cpp \
    TerminationCriterionTest.cpp \
    TopologyInfoTest.cpp \
    Vector3DTest.cpp \
    VertexCullingRegressionTest.cpp \
    VtkTest.cpp \
    DistanceFromTargetTest.cpp \
    TargetCalculatorTest.cpp 

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
        $(CPPUNIT_LNK) $(CONFIG_LDFLAGS) \
        $(TSTT_MESH_LDFLAGS) 

.SUFFIXES : .cpp .o

.cpp.o:
	@echo "Compiling $<"
	$(PREFIX)$(CXX) -c $(CXXFLAGS) $(CONFIG_CFLAGS) \
        $(MESQUITE_INCLUDE) $(CPPUNIT_INC) \
        $(TSTT_BASE_INC) $(TSTT_MESH_INC) $(TSTT_MESH_CFLAGS) \
        -c -o $@ $<
