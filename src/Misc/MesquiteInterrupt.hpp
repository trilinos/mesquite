#ifndef MESQUITEINTERRUPT_HPP
#define MESQUITEINTERRUPT_HPP

#include "Mesquite.hpp"

namespace Mesquite
{
   bool interruptFlag = false;
   
   inline bool interrupt_was_signaled()
   { return Mesquite::interruptFlag; }
   
   inline void clear_interrupt()
   { interruptFlag = false; }
}

#endif
