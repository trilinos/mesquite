#include <cppunit/extensions/TestFactoryRegistry.h>
#include "MesquiteTestRunner.hpp"

int main()
{
    // Create a test runner
  Mesquite::TestRunner runner;
  
    // Get the test suite we want to run
  CppUnit::TestFactoryRegistry &registry =
    CppUnit::TestFactoryRegistry::getRegistry("Misc");
  runner.add_test( registry.makeTest() );

    // Run the tests
  bool wasSucessful = runner.run("Darryl's Test Run");
  
    // Return 0 if there were no errors
  return wasSucessful ? 0 : 1;
}

