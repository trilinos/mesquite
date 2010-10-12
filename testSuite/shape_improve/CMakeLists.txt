#MSQ_ADD_TEST( shape_improve_test_old shape_improve_test_old_main main.cpp )
#SET_TARGET_PROPERTIES( shape_improve_test_old_main, PROPERTIES COMPILE_FLAGS -DTEST_OLD_WRAPPER )
MSQ_ADD_TEST( shape_improve_test shape_improve_test_main main.cpp )
