#include <cppunit/extensions/TestFactoryRegistry.h>
#include "MesquiteTestRunner.hpp"

int main(int argc, char **argv)
{
    // Create a test runner
  Mesquite::TestRunner runner;

    // If the user requested a specific test...
  if (argc > 1)
  {
    while (argc > 1)
    {
      argc--;
      CppUnit::TestFactoryRegistry &registry =
        CppUnit::TestFactoryRegistry::getRegistry(argv[argc]);
      runner.add_test( registry.makeTest() );
    }
  }
  else
  {
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
    
    CppUnit::TestFactoryRegistry &registry4 =
      CppUnit::TestFactoryRegistry::getRegistry("PatchDataTest");
    runner.add_test( registry4.makeTest() );

    CppUnit::TestFactoryRegistry &registry5 =
      CppUnit::TestFactoryRegistry::getRegistry("ObjectiveFunctionTest");
    runner.add_test( registry5.makeTest() );
  }
  
    // Run the tests
  bool wasSucessful = runner.run("Test Run");
  
    // Return 0 if there were no errors
  return wasSucessful ? 0 : 1;
}

