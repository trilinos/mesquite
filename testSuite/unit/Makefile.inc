UNIT_DIR = $(testdir)/unit

MSQ_TEST = $(UNIT_DIR)/msq_test

all_tests: $(MSQ_TEST)

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

TEST_OBJS_P = $(TEST_OBJS:%=$(UNIT_DIR)/%)

$(MSQ_TEST) : $(TEST_OBJS_P) $(locallibdir)/libmesquite.a
	@echo "Linking $@"
	$(PREFIX)$(LINKER) $(LDFLAGS) \
        $(TEST_OBJS_P) \
        -o $@ \
        -L$(locallibdir) -lmesquite \
        -L$(CPPUNIT_DIR)/lib -lcppunit -ldl
