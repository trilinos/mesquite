#include <cppunit/extensions/TestFactoryRegistry.h>
#include "MesquiteTestRunner.hpp"

int main()
{
    // Create a test runner
  Mesquite::TestRunner runner;
  
    // Get the test suites we want to run
  CppUnit::TestFactoryRegistry &registry =
    CppUnit::TestFactoryRegistry::getRegistry("Misc");
  runner.add_test( registry.makeTest() );

  CppUnit::TestFactoryRegistry &registry2 =
    CppUnit::TestFactoryRegistry::getRegistry("MsqMeshEntityTest");
  runner.add_test( registry2.makeTest() );

  CppUnit::TestFactoryRegistry &registry3 =
    CppUnit::TestFactoryRegistry::getRegistry("InstructionQueueTest");
  runner.add_test( registry3.makeTest() );

    // Run the tests
  bool wasSucessful = runner.run("Darryl's Test Run");
  
    // Return 0 if there were no errors
  return wasSucessful ? 0 : 1;
}

