#ifndef MESQUITE_TESTRUNNER_HPP
#define MESQUITE_TESTRUNNER_HPP

#include <cppunit/TestListener.h>
#include <string>
#include <stack>
#include <vector>
#include <iostream>

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
    virtual void startSuite(CppUnit::Test *test);
    virtual void startTest(CppUnit::Test *test);
    virtual void addFailure(const CppUnit::TestFailure &failure);
    virtual void endTest(CppUnit::Test *test);
    virtual void endSuite(CppUnit::Test *test);

      // Timer functions
    inline void push_timer(Mesquite::Timer* timer);
    inline Mesquite::Timer* pop_timer();
    
  private:
    std::vector<CppUnit::Test*> mTests;
    std::stack<Mesquite::Timer*> mTimers;
    std::ostream* mOut;
    unsigned int indentLevel;
    unsigned int numSuccesses;
    unsigned int numFailures;
    unsigned int numExceptions;
    static const unsigned char INDENT_SIZE;
  };
} // namespace Mesquite

#endif
