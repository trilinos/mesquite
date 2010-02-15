if ENABLE_IMESH
  IMESH_TEST = imesh
else
  IMESH_TEST =
endif

if ENABLE_IGEOM
  IGEOM_TEST = igeom
else
  IGEOM_TEST =
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
	  idft_time \
	  headers \
	  2d_target \
	  2d_metrics\
	  2d_formulation \
	  synchronous \
          high_aspect_ratio \
	  higher_order \
	  viscous_cfd \
	  size_adapt_shape \
	  $(IMESH_TEST) $(IGEOM_TEST) \
          $(EXO_TEST) $(UNIT_TEST)

EXTRA_DIST = CMakeLists.txt

# Override all-recursive so that "make" with no target 
# skips subdirectories in testSuite.
all-recursive:
	@:
