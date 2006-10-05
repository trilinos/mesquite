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
#ifndef MESQUITE_TESTRUNNER_HPP
#define MESQUITE_TESTRUNNER_HPP

#include <cppunit/TestListener.h>
#include <cppunit/TestSuite.h>
#include "Mesquite.hpp"
#ifdef MSQ_USE_OLD_STD_HEADERS
#include <string.h>
#include <stack.h>
#include <vector.h>
#else
#include <string>
#include <stack>
#include <vector>
#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else
#include <iostream>
#endif

namespace Mesquite
{
  class Timer;

/*!
 * \brief A class the runs cppunit tests, outputs results in an organized manner.
 *
 * The test runner manages the life cycle of the added tests.
 *
 * Here is an example of use:
 *
 * \code
 * Mesquite::TestRunner runner;
 * runner.addTest( ExampleTestCase::suite() );
 * runner.run( "Darryl's Test Run" );    // Run all tests and wait
 * \endcode
 *
 */
  class CPPUNIT_API TestRunner : protected CppUnit::TestListener
  {
  public:
    TestRunner();
    virtual ~TestRunner();
    
    void add_test(CppUnit::Test *test);
    virtual bool run(const std::string& name_of_run,
                     std::ostream& out_stream = std::cout);
  protected:
    void delete_all_tests();
    const std::string running_test_prefix();
    inline void indent();
    
      // TestListener functions
    virtual void startSuite(CppUnit::TestSuite *suite);
    virtual void startTest(CppUnit::Test *test);
    virtual void addFailure(const CppUnit::TestFailure &failure);
    virtual void endTest(CppUnit::Test *test);
    virtual void endSuite(CppUnit::TestSuite *suite);

      // Timer functions
    inline void push_timer(Mesquite::Timer* timer);
    inline Mesquite::Timer* pop_timer();
    
  private:
    msq_std::vector<CppUnit::Test*> mTests;
    msq_std::stack<Mesquite::Timer*> mTimers;
    msq_std::stack<int> failureCounters;
    msq_std::vector<CppUnit::TestSuite::Key> completedSuites;
    msq_std::vector<msq_std::string> failedTestNames;
    msq_stdio::ostream* mOut;
    CppUnit::TestResult* myResult;
    unsigned int indentLevel;
    unsigned int numSuccesses;
    unsigned int numFailures;
    unsigned int numExceptions;
    static const unsigned char INDENT_SIZE;
  };
} // namespace Mesquite

#endif
