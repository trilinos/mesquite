#ifndef MESQUITEINTERRUPT_HPP
#define MESQUITEINTERRUPT_HPP

#include "Mesquite.hpp"
#ifdef ENABLE_INTERRUPT
#include <signal.h>
#include "MsqMessage.hpp"
#endif

namespace Mesquite
{
  class MesquiteInterrupt
  {
  public:
#ifdef ENABLE_INTERRUPT
    static void catch_interrupt( bool yesno );
#endif
    static inline bool interrupt_was_signaled()
      { return interruptFlag; }
   
    static inline void clear_interrupt()
      { interruptFlag = false; }
     //michael should this be 'volatile'
    static bool interruptFlag;
  };
  
  
}

#endif
