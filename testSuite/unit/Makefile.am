@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/libmesquite.la
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@


if WITH_CPPUNIT
  check_PROGRAMS = msq_test
  msq_test_DEPENDENCIES = $(top_builddir)/libmesquite.la
  TESTS = $(check_PROGRAMS)
  CPPUNIT_CPPFLAGS_SET = @CPPUNIT_CPPFLAGS@
  CPPUNIT_LDFLAGS_SET = @CPPUNIT_LDFLAGS@
else
  TESTS = no_cpp_unit
  CPPUNIT_CPPFLAGS_SET = 
  CPPUNIT_LDFLAGS_SET = 
endif
AM_CPPFLAGS += $(CPPUNIT_CPPFLAGS_SET)
AM_LDFLAGS += $(CPPUNIT_LDFLAGS_SET)

no_cpp_unit:
	@echo "echo '###################################################'"  > $@
	@echo "echo '#### Mesquite configured w/out CppUnit support ####'" >> $@
	@echo "echo '####          Cannot build unit tests          ####'" >> $@
	@echo "echo '###################################################'" >> $@
	@echo "exit 0" >> $@
	@chmod a+x $@

if WITH_TEST_IMESH
  @MESQUITE_IMESH_MAKE_INCLUDE@
  IMESH_TEST_SRC = iMeshTest.cpp
  IMESH_LDADD = $(IMESH_LIBS)
else
  IMESH_TEST_SRC =
  IMESH_LDADD = 
endif
LDADD += $(IMESH_LDADD)

if WITH_EXODUS
  EXODUS_TEST_SRC = ExodusTest.cpp
else
  EXODUS_TEST_SRC = 
endif
  
msq_test_SOURCES = msq_test_main.cpp \
	MesquiteTestRunner.hpp \
	ObjectiveFunctionTests.hpp \
	PatchDataInstances.hpp \
	QualityMetricTester.hpp \
	QualityMetricTester.cpp \
	UnitUtil.hpp \
	ArrayMeshTest.cpp \
	AspectRatioGammaTest.cpp \
	AveragingQMTest.cpp \
	BoundedCylinderDomainTest.cpp \
	CachingTargetTest.cpp \
	CompositeMetricTest.cpp \
	CompositeOFTest.cpp \
	ConditionNumberTest.cpp \
	CylinderDomainTest.cpp \
        DomainClassifierTest.cpp \
	EdgeLengthMetricTest.cpp \
	ElemSampleQMTest.cpp \
	ExtraDataTest.cpp \
	FileTokenizerTest.cpp \
        GeomPrimTest.cpp \
	IdealElementTest.cpp \
	IdealTargetTest.cpp \
	IdealWeightMeanRatioTest.cpp \
	IdealWeightInverseMeanRatioTest.cpp \
	InstructionQueueTest.cpp \
	I_DFT_Test.cpp \
	TMPQualityMetricTest.cpp \
	LinearMappingFunctionTest.cpp \
	LVQDTargetTest.cpp \
	Matrix3DTest.cpp \
	MeshInterfaceTest.cpp \
	MeshImplTest.cpp \
	MsqFreeVertexIndexIteratorTest.cpp \
	MsqHessianTest.cpp \
	MsqMatrixTest.cpp \
	MsqMeshEntityTest.cpp \
	MsqVertexTest.cpp \
	NumericalOFTest.cpp \
	ObjectiveFunctionTest.cpp \
	PatchDataTest.cpp \
	PatchDataTestNormals.cpp \
	PatchSetTest.cpp \
	PlanarGeometryTest.cpp \
	PMeanPMetricTest.cpp \
	PMeanPTemplateTest.cpp \
	QuadLagrangeShapeTest.cpp \
	QualityAssessorTest.cpp \
	QualityMetricTest.cpp \
	SamplePointsTest.cpp \
	SlaveBoundaryVerticesTest.cpp \
	SphericalGeometryTest.cpp \
	StdDevTemplateTest.cpp \
	SymMatrix3DTest.cpp \
	TagVertexMeshTest.cpp \
	Target2DTest.cpp \
	Target3DTest.cpp \
	TargetReadWriteTest.cpp \
	TerminationCriterionTest.cpp \
	TetLagrangeShapeTest.cpp \
	TriLagrangeShapeTest.cpp \
	TopologyInfoTest.cpp \
	UntangleBetaTest.cpp \
	Vector3DTest.cpp \
	VertexConditionNumberTest.cpp \
	VertexCullingRegressionTest.cpp \
	VtkTest.cpp \
	XYRectangleTest.cpp \
	BCDTest.cpp \
	Target2DSurfOrientTest.cpp \
        $(IMESH_TEST_SRC) \
        $(EXODUS_TEST_SRC)

CLEANFILES = no_cpp_unit
MOSTLYCLEANFILES = hexagon.vtk \
		   LPtoP-hex-bcd.vtk \
		   LPtoP-hex-gbl.vtk \
		   LPtoP-tet-bcd.vtk \
		   LPtoP-tet-gbl.vtk \
		   PMeanP-hex-bcd.vtk \
		   PMeanP-hex-gbl.vtk \
		   PMeanP-tet-bcd.vtk \
		   PMeanP-tet-gbl.vtk \
		   domain_ideal.vtk \
		   domain_ref.vtk \
		   none_ideal.vtk \
		   none_ref.vtk \
		   target_ideal.vtk \
		   target_ref.vtk

		   
