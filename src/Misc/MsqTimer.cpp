#include "MsqTimer.hpp"

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

double Mesquite::Timer::since_birth()
{
  return now() - atBirth;
}
