#include "MsqTimer.hpp"

// Create the global collection of stop watches
Mesquite::StopWatchCollection Mesquite::GlobalStopWatches;

// Unfortunately, many Unix platforms don't have good support
// for the ANSI standard clock() function.  That means we have
// to implement this class differently for WIN32 and Unix.  If
// your platform wants to use clock() to implement this class,
// define USE_CLOCK_TIMER in your compile line for this file.
// If USE_CLOCK_TIMER is not defined, the class is implemented
// using times() instead.
#ifdef USE_CLOCK_TIMER
#ifdef USE_C_PREFIX_INCLUDES
#include <ctime>
#else
#include <time.h>	
#endif
static inline double now()
{ return (double)clock() / (double)CLOCKS_PER_SEC; }

#else // If we aren't using clock()

#include <sys/param.h>
#include <sys/times.h>
#include <time.h>	
#ifdef SOLARIS
#include <unistd.h>
#endif

#ifndef HZ
#ifdef CLK_TCK
#define HZ CLK_TCK
#else
#ifdef CLOCKS_PER_SEC
#define HZ CLOCKS_PER_SEC
#else
#define HZ 60
#endif
#endif
#endif

static inline double now()
{
  tms current;
  times( &current );
  return ((double)(current.tms_utime  +
                   current.tms_stime  +
                   current.tms_cutime +
                   current.tms_cstime)) / (double)HZ;
}
#endif // End of times() based implementation

Mesquite::Timer::Timer() 
    : atBirth(now())
{
  atLastCheck = atBirth;
}

double Mesquite::Timer::since_last_check()
{
  double right_now = now();
  double rv = right_now - atLastCheck;
  atLastCheck = right_now;
  return rv;
}

double Mesquite::Timer::since_birth() const
{
  return now() - atBirth;
}

void Mesquite::StopWatch::start()
{
  if (!isRunning)
  {
    isRunning = true;
    timeAtLastStart=now();
    ++numStarts;
  }
}

void Mesquite::StopWatch::stop()
{
  if (isRunning)
  {
    isRunning = false;
    totalTime += now() - timeAtLastStart;
  }
}

void Mesquite::StopWatch::reset()
{
  isRunning=false;
  totalTime=0;
  numStarts=0;
}

double Mesquite::StopWatch::total_time() const
{
  double rv = totalTime;
  if (isRunning)
    rv += now() - timeAtLastStart;
  return rv;
}

Mesquite::StopWatchCollection::Key Mesquite::StopWatchCollection::add(
  const std::string &name,
  bool fail_if_exists)
{
    // Don't allow empty name
  if (name == "")
    return 0;
  
  Key key = get_key(name);
  
    // If the named stopwatch doesn't exist...
  if (!key)
  {
      // See if there is an unused existing stopwatch
    int i;
    for (i = 0; i < mEntries.size(); i++)
    {
      if (mEntries[i].first == "")
      {
        mEntries[i].first = name;
        mEntries[i].second.reset();
        break;
      }
    }
      // If not, create a new one
    if (i == mEntries.size())
    {
      mEntries.push_back(std::pair<std::string, StopWatch>(name, StopWatch()));
    }
    key = i+1;
  }
    // If it already existed...
  else if (fail_if_exists)
    key = 0;
  
  return key;
}


Mesquite::StopWatchCollection::Key Mesquite::StopWatchCollection::get_key(
  const std::string &name) const
{
  Key key = 0;
  
  for (int i = 0; i < mEntries.size(); i++)
  {
    if (mEntries[i].first == name)
    {
      key = i + 1;
      break;
    }
  }

  return key;
}

void Mesquite::StopWatchCollection::remove(
  const Mesquite::StopWatchCollection::Key key)
{
    // Get rid of anything at the end of the list
  if (key == mEntries.size())
  {
    mEntries.pop_back();
    while (!mEntries.empty() && mEntries.back().first == "")
    {
      mEntries.pop_back();
    }
  }
  
  else if (key > 0 && key < mEntries.size())
  {
      // If in the middle of the list, set its name to ""
    mEntries[key-1].first = "";
  }
}


void Mesquite::StopWatchCollection::start(
  const Mesquite::StopWatchCollection::Key key)
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    mEntries[key-1].second.start();
}

void Mesquite::StopWatchCollection::stop(
  const Mesquite::StopWatchCollection::Key key)
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    mEntries[key-1].second.stop();
}

void Mesquite::StopWatchCollection::reset(
  const Mesquite::StopWatchCollection::Key key)
{
  if (key > 0 &&
      key <= mEntries.size())
    mEntries[key-1].second.reset();
}


double Mesquite::StopWatchCollection::total_time(
  const Mesquite::StopWatchCollection::Key key) const
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    return mEntries[key-1].second.total_time();
  else
    return 0.0;
}

int Mesquite::StopWatchCollection::number_of_starts(
  const Mesquite::StopWatchCollection::Key key) const
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    return mEntries[key-1].second.number_of_starts();
  else
    return 0;
}
/*! Fills a vector of StopWatchCollection::Key in which the Keys are ordered
  by the associated StopWatch's total_time.  The key associated with the
  largest total_time StopWatch is in the first position of the vector.  The
  key associated with the smallest total_time StopWatch is in the last
  position of the vector.*/
void Mesquite::StopWatchCollection::get_keys_sorted_by_time(
  std::vector<Key> &sorted_keys)
{
  int num_watches=mEntries.size();
  int *sorted_indices=new int[num_watches];
  int i=0;
  int counter=0;
  for(i=0;i<num_watches;++i){
    sorted_indices[i]=0;
  }
  double current_max;
  int index_to_max;
    //While we haven't added all of the Keys to the vector
  while(counter<num_watches){
    current_max=-1;
    index_to_max=-1;
      //loop over the times and find the largest remaining
    for(i=0;i<num_watches;++i){
      if(mEntries[i].second.total_time()>current_max && sorted_indices[i]==0){
        current_max=mEntries[i].second.total_time();
        index_to_max=i;
      }
    }
      //Add the key associated with index_to_max and any subsequent
      //keys which are associated with a StopWatch that has a total
      //time equal to current_max;
    for(i=index_to_max;i<=num_watches;++i){
      if(mEntries[i].second.total_time()>=current_max && sorted_indices[i]==0)
      {
        counter++;
        sorted_indices[i]=counter;
        sorted_keys.push_back(i+1);
      }
    }
  }
    //clean up
  delete[] sorted_indices;
}

      
  

