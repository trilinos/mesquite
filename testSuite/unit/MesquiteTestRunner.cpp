#include "MesquiteTestRunner.hpp"
#include "MsqTimer.hpp"
#include "cppunit/Test.h"
#include "cppunit/TestResult.h"
#include "cppunit/TestFailure.h"
#include "cppunit/Exception.h"

const unsigned char Mesquite::TestRunner::INDENT_SIZE = 2;

void Mesquite::TestRunner::indent()
{
  unsigned int local_indent = indentLevel;
  while (local_indent--)
    *mOut << ' ';
}

void Mesquite::TestRunner::push_timer(Mesquite::Timer* timer)
{
  mTimers.push(timer);
}

Mesquite::Timer* Mesquite::TestRunner::pop_timer()
{
  Mesquite::Timer* rv = mTimers.top();
  mTimers.pop();
  return rv;
}

Mesquite::TestRunner::TestRunner() 
    : indentLevel(0)
{}

Mesquite::TestRunner::~TestRunner()
{
  delete_all_tests();
}

void Mesquite::TestRunner::delete_all_tests()
{
  for (std::vector<CppUnit::Test*>::iterator iter = mTests.begin();
       iter != mTests.end();
       ++iter)
  {
    delete *iter;
  }
  mTests.clear();
}

bool Mesquite::TestRunner::run(const std::string& name_of_run,
                               std::ostream& out_stream)
{
  mOut = &out_stream;
  *mOut << "Start of test run: " << name_of_run << "\n\n";
  
    // Run each test
  CppUnit::TestResult result;
  result.addListener(this);
  for (std::vector<CppUnit::Test*>::iterator iter = mTests.begin();
       iter != mTests.end();
       ++iter)
    (*iter)->run(&result);
  
  *mOut << "\nEnd of test run: " << name_of_run << std::endl;

  return true;
}

void Mesquite::TestRunner::add_test(CppUnit::Test *test)
{
  mTests.push_back(test);
}

// This function is called by the TestRunner just before
// a test begins.
void Mesquite::TestRunner::startTest(CppUnit::Test *test)
{
    // Indent
  indent();
  
    // Output a header
  *mOut << "Beginning of Test : " << test->getName() << std::endl;
  
    // increase the indent level
  indentLevel += TestRunner::INDENT_SIZE;

    // Start a timer
  push_timer(new Mesquite::Timer);
}

// This function is called by the TestRunner just before
// a suite begins.
void Mesquite::TestRunner::startSuite(CppUnit::Test *test)
{
    // Indent
  indent();
  
    // Output a header
  *mOut << "Beginning of Test Suite : " << test->getName() << std::endl;
  
    // increase the indent level
  indentLevel += TestRunner::INDENT_SIZE;

    // Add a timer
  push_timer(new Mesquite::Timer);
}

// This function is called if a test fails, either
// intentionally or unintentionally.
void Mesquite::TestRunner::addFailure(const CppUnit::TestFailure &failure)
{
    // Indent
  indent();
  
  if (failure.isError())
  {
    *mOut << "***ERROR***\n";
    indent();
    if (failure.thrownException())
    {
      *mOut << "Unexpected exception : " << failure.thrownException()->what()
            << std::endl;
    }
    else
    {
      *mOut << "Unexpected test failure" << std::endl;
    }
  }
  else
  {
    *mOut << "Failed assertion : " << std::endl;
  }
  
}

// This function is called just after a test completes.
void Mesquite::TestRunner::endTest(CppUnit::Test *test)
{
    // Pop the timer
  Mesquite::Timer *timer = pop_timer();
  double elapsed_time = timer->since_birth();
  delete timer;
  
    // Decrease the indent level
  indentLevel -= TestRunner::INDENT_SIZE;
  
    // Output a footer
  indent();
  *mOut << "Elapsed time: " << elapsed_time << " seconds\n";
  indent();
  *mOut << "End of Test : " << test->getName() << std::endl;
}

// This function is called just after a test completes.
void Mesquite::TestRunner::endSuite(CppUnit::Test *test)
{
    // Pop the timer
  Mesquite::Timer *timer = pop_timer();
  double elapsed_time = timer->since_birth();
  delete timer;
  
    // Decrease the indent level
  indentLevel -= TestRunner::INDENT_SIZE;
  
    // Output a footer
  indent();
  *mOut << "Elapsed time: " << elapsed_time << " seconds\n";
  indent();
  *mOut << "End of Test Suite : " << test->getName() << std::endl;
}
