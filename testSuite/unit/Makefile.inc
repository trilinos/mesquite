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

TEST_SRCS = \
    BoundedCylinderDomainTest.cpp \
    CylinderDomainTest.cpp \
    ExodusTest.cpp \
    FileTokenizerTest.cpp \
    InstructionQueueTest.cpp \
    Matrix3DTest.cpp \
    MeshInterfaceTest.cpp \
    MsqFreeVertexIndexIteratorTest.cpp \
    MsqHessianTest.cpp \
    MsqMeshEntityTest.cpp \
    MsqVertexTest.cpp \
    ObjectiveFunctionTest.cpp \
    PatchDataTest.cpp \
    PatchDataTestNormals.cpp \
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

$(MSQ_TEST) : $(TEST_OBJS_P) testSuite/unit/msq_test_main.o $(locallibdir)/libmesquite.a
	@echo "Linking $@"
	$(PREFIX)$(LINKER) $(LDFLAGS) \
	$(TEST_OBJS_P) testSuite/unit/msq_test_main.o \
	-o $@ \
	-L$(locallibdir) -lmesquite \
	$(CPPUNIT_LNK) $(CONFIG_LDFLAGS) \
        $(TSTTM_SERVER_LIBS) $(TSTTB_SERVER_LIBS)

testSuite/unit/qt_test : $(TEST_OBJS_P) testSuite/unit/qt_test_main.o $(locallibdir)/libmesquite.a
	@echo "Linking $@"
	$(PREFIX)$(LINKER) $(LDFLAGS) \
	$(TEST_OBJS_P) testSuite/unit/qt_test_main.o \
	-o $@ \
	-L$(locallibdir) -lmesquite \
	$(CPPUNIT_LNK) $(CONFIG_LDFLAGS) \
	-lqttestrunner -L$(QTDIR)/lib -lqt-mt \
        $(TSTTM_SERVER_LIBS) $(TSTTB_SERVER_LIBS)

.SUFFIXES : .cpp .o

qt_test_main.o: qt_test_main.cpp
	$(PREFIX) $(CXX) -c $(CXXFLAGS) $(CONFIG_CFLAGS) \
        $(MESQUITE_INCLUDE) $(CPPUNIT_INC) -I$(QTDIR)/include \
        $(TSTT_BASE_INC) $(TSTT_MESH_INC) \
        -c -o $@ $<

.cpp.o:
	@echo "Compiling $<"
	$(PREFIX)$(CXX) -c $(CXXFLAGS) $(CONFIG_CFLAGS) \
        $(MESQUITE_INCLUDE) $(CPPUNIT_INC) \
        $(TSTT_BASE_INC) $(TSTT_MESH_INC) \
        -c -o $@ $<
