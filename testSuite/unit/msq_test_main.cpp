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
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/Outputter.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestFailure.h>
#include <assert.h>
#include <stdio.h>

class CPPUNIT_API SummaryOutput : public CppUnit::Outputter 
{
  public:
    SummaryOutput( FILE* file, CppUnit::TestResultCollector* result ) 
      : file_(file), results_(result) {}
    void write();
  private:
    FILE* file_;
    CppUnit::TestResultCollector* results_;
};


int main(int argc, char **argv)
{
  CppUnit::TextUi::TestRunner runner;
  int firsttest = 1;

  if (argc > 1 && !strcmp(argv[1],"-s"))
  {
    FILE* file = fopen(argv[2],"w");
    if (!file) 
    {
      perror( argv[2] );
      exit(1);
    }
    runner.setOutputter( new SummaryOutput( file, &runner.result() ) );
    firsttest += 2;
  }

    // If the user requested a specific test...
  if (argc > firsttest)
  {
    while (argc > firsttest)
    {
      argc--;
      CppUnit::TestFactoryRegistry &registry =
        CppUnit::TestFactoryRegistry::getRegistry(argv[argc]);
      runner.addTest( registry.makeTest() );
    }
    
  }
  else
  {
     CppUnit::Test* test;
     test = CppUnit::TestFactoryRegistry::getRegistry("Unit").makeTest();
     runner.addTest( test );
     test = CppUnit::TestFactoryRegistry::getRegistry("Regression").makeTest();
     runner.addTest( test );
  }
  
    // Return 0 if there were no errors
  return !runner.run( );
}

void SummaryOutput::write()
{
  CppUnit::TestResultCollector::TestFailures fails = results_->failures();
  CppUnit::TestResultCollector::Tests tests = results_->tests();
  
  CppUnit::TestResultCollector::TestFailures::const_iterator f_iter = fails.begin();
  CppUnit::TestResultCollector::Tests::const_iterator t_iter;
  
  fprintf(file_,"****Tests Run:\n");
  for (t_iter = tests.begin(); t_iter != tests.end(); ++t_iter)
    fprintf(file_, "%s\n", (*t_iter)->getName().c_str());

  fprintf(file_,"****Tests Failed:\n");
  for (f_iter = fails.begin(); f_iter != fails.end(); ++f_iter)
    fprintf(file_, "%s\n", (*f_iter)->failedTestName().c_str());
 
}

