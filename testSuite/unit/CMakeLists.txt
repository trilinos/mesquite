SET( unit_tests
     ArrayMeshTest 
     AspectRatioGammaTest 
     AveragingQMTest
     BoundedCylinderDomainTest
     CachingTargetTest
     CompositeMetricTest
     CompositeOFTest
     ConditionNumberTest
     CylinderDomainTest
     DomainClassifierTest
     EdgeLengthMetricTest
     EdgeLengthQualityMetricTest
     ElemSampleQMTest
     ExtraDataTest
     FileTokenizerTest
     GeomPrimTest
     HigherOrderTest
     IdealElementTest
     IdealTargetTest
     IdealWeightMeanRatioTest
     IdealWeightInverseMeanRatioTest
     InstructionQueueTest
     TMPQualityMetricTest
     LinearMappingFunctionTest
     LVQDTargetTest
     MappingFunctionTest
     Matrix3DTest
     MeshInterfaceTest
     MeshImplTest
     MsqFreeVertexIndexIteratorTest
     MsqHessianTest
     MsqMatrixTest
     MsqMeshEntityTest
     MsqVertexTest
     NodeSetTest
     NumericalOFTest
     ObjectiveFunctionTest
     PatchDataTest
     PatchDataTestNormals
     PatchSetTest
     PlanarGeometryTest
     PMeanPMetricTest
     PMeanPTemplateTest
     QuadLagrangeShapeTest
     QualityAssessorTest
     QualityMetricTest
     SlaveBoundaryVerticesTest
     SphericalGeometryTest
     StdDevTemplateTest
     SymMatrix3DTest
     TagVertexMeshTest
     Target2DSurfOrientTest
     Target2DTest
     Target3DTest
     TargetReadWriteTest
     TerminationCriterionTest
     TetLagrangeShapeTest
     TriLagrangeShapeTest
     TopologyInfoTest
     TMPDerivsTest
     UntangleBetaTest
     Vector3DTest
     VertexConditionNumberTest
     VertexCullingRegressionTest
     VtkTest
     XYRectangleTest
     BCDTest )

SET( unit_libs ${CPPUNIT_LIBRARY} )

IF (${PACKAGE_NAME}_ENABLE_IMESH)
  SET(unit_tests ${unit_tests} iMeshTest)
  SET(unit_libs ${unit_libs} ${IMESH_LIBRARY})
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_ExodusII)
  SET(unit_tests ${unit_tests} ExodusTest)
  SET(unit_libs ${unit_libs} ${EXODUS_LIBRARY})
ENDIF()

SET( unit_srcs
     msq_test_main.cpp 
     QualityMetricTester.cpp 
     MesquiteTestRunner.hpp 
     ObjectiveFunctionTests.hpp 
     PatchDataInstances.hpp 
     QualityMetricTester.hpp 
     UnitUtil.hpp )

FOREACH( test ${unit_tests} )
  SET( unit_srcs ${unit_srcs} ${test}.cpp )
ENDFOREACH()

IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
  PACKAGE_ADD_EXECUTABLE( msq_test NOEXEPREFIX NOEXESUFFIX SOURCES ${unit_srcs} )
  TARGET_LINK_LIBRARIES( msq_test mesquite ${CPPUNIT_LIBRARY} )
  FOREACH( test ${unit_tests} )
    PACKAGE_ADD_TEST( msq_test NAME ${test} ARGS ${test} )
  ENDFOREACH()
ELSE()
  IF (${PACKAGE_NAME}_ENABLE_TESTS)
    ADD_EXECUTABLE( msq_test ${unit_srcs} )
    TARGET_LINK_LIBRARIES( msq_test mesquite ${CPPUNIT_LIBRARY} )
    FOREACH( test ${unit_tests} )
      ADD_TEST( ${test} ${EXECUTABLE_OUTPUT_PATH}/msq_test ${test} )
    ENDFOREACH()
  ENDIF()
ENDIF()

INCLUDE_DIRECTORIES( ${CPPUNIT_INCLUDES} )
SET_SOURCE_FILES_PROPERTIES(${unit_srcs} COMPILE_FLAGS "${TEST_DEFINES} '-DSRCDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"'")
