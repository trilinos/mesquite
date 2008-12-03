if WITH_TEST_IMESH
  IMESH_TEST = imesh
else
  IMESH_TEST =
endif

if WITH_EXODUS
  EXO_TEST = convert 
else
  EXO_TEST = 
endif

if WITH_CPPUNIT
  UNIT_TEST = unit
else
  UNIT_TEST =
endif

SUBDIRS = test_1 \
	  ActiveSetTest \
	  laplacian_test\
	  laplacian_wrapper_test\
	  untangle_test\
	  analytical_grad_3D\
	  algorithm_test\
	  simple_hybrid_test\
	  wrapper_tests\
	  tutorial\
	  transform \
	  pyramid \
	  wedge \
	  Guides843 \
	  deforming \
	  rezone \
	  i_dft_with_culling \
	  idft_time \
	  test_DFT \
	  headers \
	  2d_target \
	  2d_metrics\
	  2d_formulation \
	  synchronous \
          high_aspect_ratio \
	  $(IMESH_TEST) $(EXO_TEST) $(UNIT_TEST)

# Override all-recursive so that "make" with no target 
# skips subdirectories in testSuite.
all-recursive:
	@:
