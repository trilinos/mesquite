#ifndef MESQUITE_TIMER_HPP
#define MESQUITE_TIMER_HPP

#include <vector>
#include <utility>
#include <string>

namespace Mesquite
{
  class Timer
  {
  public:
    Timer();
    
    double since_last_check(); //- return time in seconds since last
                               //- call to since_last_check().  Note that
                               //- calling since_birth() doesn't count
                               //- as a check.  The first time this function
                               //- is called, it returns the time since birth.
    
    double since_birth() const;//- return time in seconds since
                               //- object was created.
    
  private:
    const double atBirth;      //- Time at birth
    double atLastCheck;        //- Time at last call to since_last_check()
  };
  
  
  class StopWatch
  {
  public:
      // Creates the stopwatch.  The stopwatch is stopped
      // until start() is called.
    StopWatch() 
        :isRunning(false), totalTime(0)
      {}
    
      // Starts the stopwatch.  If it was already running,
      // this function does nothing.
    void start();
    
      // Stops the stopwatch.  If it was not already running,
      // this function does nothing
    void stop();
    
      // Stops the stopwatch and resets the total_time() to zero.
    void reset();
    
      // Returns the total accumulated time.  If the stopwatch
      // is currently running, the time between the last start()
      // and the current time IS included in total_time().
    double total_time() const;
    
  private:
    bool isRunning;
    double timeAtLastStart;
    double totalTime;
  };
  
  class StopWatchCollection
  {
  public:
    typedef size_t Key;
    
      // Create a new collection
    StopWatchCollection()
      {}
    
      // Add a stopwatch to the collection.  Returns a non-zero
      // StopWatchCollection::Key if it succeeds, zero if it fails.
      // If a StopWatch with the given name already exists in the
      // collection, the Key of the existing StopWatch is returned
      // if 'fail_if_exists' is false, or zero is returned if
      // 'fail_if_exists' is true.
    Key add(const std::string &name, bool fail_if_exists = true);
    
      // Gets the Key for an existing stopwatch.  If a stopwatch
      // with the given name does not exist, function returns zero.
    Key get_key(const std::string &name) const;
    
      // Remove a specific stopwatch.
    void remove(const Key key);
    void remove(const std::string &name)
      { remove(get_key(name)); }
    
      // start a specific stopwatch
    void start(const Key key);
    void start(const std::string &name)
      { start(get_key(name)); }
    
      // stop a specific stopwatch
    void stop(const Key key);
    void stop(const std::string &name)
      { stop(get_key(name)); }
    
      // reset a specific stopwatch
    void reset(const Key key);
    void reset(const std::string &name)
      { reset(get_key(name)); }
    
      // Get the total time for a specific stopwatch, zero if
      // the stopwatch doesn't exist.
    double total_time(const Key key) const;
    double total_time(const std::string &name) const
      { return total_time(get_key(name)); }
    
  private:
    std::vector< std::pair<std::string, StopWatch> > mEntries;
  };
  
    // A stopWatchCollection available anywhere
  extern Mesquite::StopWatchCollection GlobalStopWatches;
}

#endif
