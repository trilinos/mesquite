#ifndef MESQUITE_TIMER_HPP
#define MESQUITE_TIMER_HPP

namespace Mesquite
{
  class Timer
  {
  public:
    Timer();
    
    double since_last_check(); //- return CPU time in seconds since last
                               //- call to since_last_check().  Note that
                               //- calling since_birth() doesn't count
                               //- as a check.
    
    double since_birth();      //- return CPU time in seconds since
                               //- object was created.
    
  private:
    const double atBirth;      //- Time at birth
    double atLastCheck;        //- Time at last call to since_last_check()
  };
}

#endif
