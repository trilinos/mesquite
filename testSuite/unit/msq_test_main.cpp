/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/SignalException.h>
#include <signal.h>
#include "MesquiteTestRunner.hpp"
#include "MsqMessage.hpp"
int main(int argc, char **argv)
{
    // Create a test runner
  Mesquite::TestRunner runner;
  bool wasSuccessful = false;
  
  CppUnit::SignalException::throw_on_signal(SIGSEGV, false);
  
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

    CppUnit::TestFactoryRegistry &registry6 =
      CppUnit::TestFactoryRegistry::getRegistry("MeshSetTest");
    runner.add_test( registry6.makeTest() );

    CppUnit::TestFactoryRegistry &registry7 =
      CppUnit::TestFactoryRegistry::getRegistry("MsqVertexTest");
    runner.add_test( registry7.makeTest() );

     CppUnit::TestFactoryRegistry &registry8 =
      CppUnit::TestFactoryRegistry::getRegistry("MsqFreeVertexIndexIteratorTest");
    runner.add_test( registry8.makeTest() );

    CppUnit::TestFactoryRegistry &registry9 =
      CppUnit::TestFactoryRegistry::getRegistry("QualityMetricTest");
    runner.add_test( registry9.makeTest() );

    CppUnit::TestFactoryRegistry &registry10 =
      CppUnit::TestFactoryRegistry::getRegistry("AomdVtkTest");
    runner.add_test( registry10.makeTest() );
  }
  
    // Run the tests
  wasSuccessful = runner.run("Test Run");

    // Return 0 if there were no errors
  return wasSuccessful ? 0 : 1;
}

