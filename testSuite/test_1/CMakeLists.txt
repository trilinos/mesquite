
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Test_1
  SOURCES main.cpp
  COMM serial mpi
  )
