#ifndef MESQUITE_TIMER_HPP
#define MESQUITE_TIMER_HPP

#ifdef WIN32
#pragma warning ( 4 : 4786)
#endif

#include <vector>
#include <utility>
#include <string>

namespace Mesquite
{
  class Timer
  {
  public:
    Timer();

    void reset();//resets the timer as if it were just created
    
    double since_last_check(); //- return time in seconds since last
                               //- call to since_last_check().  Note that
                               //- calling since_birth() doesn't count
                               //- as a check.  The first time this function
                               //- is called, it returns the time since birth.
    
    double since_birth() const;//- return time in seconds since
                               //- object was created.
    
  private:
    double atBirth;      //- Time at birth
    double atLastCheck;        //- Time at last call to since_last_check()
  };
  
  
  class StopWatch
  {
  public:
      // Creates the stopwatch.  The stopwatch is stopped
      // until start() is called.
    StopWatch() 
        :isRunning(false), totalTime(0.0), numStarts(0)
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
      
      /*! \brief Returns the number of times this StopWatch has
        been started.*/
    int number_of_starts() const{
        return numStarts;
      }

    
  private:
    bool isRunning;
    double timeAtLastStart;
    double totalTime;
    int numStarts;
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

      //!Gets the string associated with a key
    std::string get_string(const Key key){
        return mEntries[key-1].first;}
      //!Gets the string associated with a key      
    void get_string(const Key key, std::string &new_string){
      new_string=mEntries[key-1].first;}
    
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
      // Get the number of times a StopWatch was started.
    int number_of_starts(const Key key) const;
    int number_of_starts(const std::string &name) const
      { return number_of_starts(get_key(name));}
    
      //Gets the number of stop watches in the collection
    int number_of_stop_watches(){
      return mEntries.size();}

    void get_keys_sorted_by_time(std::vector<Key> &sorted_keys);
    
    
  private:
    std::vector< std::pair<std::string, StopWatch> > mEntries;
  };
  
    // A stopWatchCollection available anywhere
  extern Mesquite::StopWatchCollection GlobalStopWatches;
}

#ifdef USE_FUNCTION_TIMERS
#define FUNCTION_TIMER_START(name) \
   static StopWatchCollection::Key local_func_abd = GlobalStopWatches.add(name); \
   GlobalStopWatches.start(local_func_abd); 

#define FUNCTION_TIMER_END() \
   GlobalStopWatches.stop(local_func_abd); 
#else
#define FUNCTION_TIMER_START(name){}
#define FUNCTION_TIMER_END() {}
#endif


#endif
