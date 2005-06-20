UNIT_DIR = $(testdir)/unit

UNIT_TEST_EXE = msq_test

MSQ_TEST = $(UNIT_DIR)/$(UNIT_TEST_EXE)

MSQ_RUN_TEST = $(UNIT_DIR)/run_target
.Phoney: $(MSQ_TARGET)

ifeq ($(CPPUNIT_LNK),)
$(MSQ_RUN_TEST):
	@echo "====================================================" 
	@echo "Reconfigure with CPPUnit support to build unit tests" 
	@echo "====================================================" 
else
$(MSQ_RUN_TEST): $(MSQ_TEST)
	cd $(UNIT_DIR) && ./$(UNIT_TEST_EXE)
all_tests: $(MSQ_TEST)
endif

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
