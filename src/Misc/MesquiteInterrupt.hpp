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
#ifndef MESQUITEINTERRUPT_HPP
#define MESQUITEINTERRUPT_HPP

#include "Mesquite.hpp"
#ifdef ENABLE_INTERRUPT
#include <signal.h>
#include "MsqMessage.hpp"
#endif

//do this because sighandler_t is not defined on all platforms
extern "C"
{
  typedef void (*msq_sig_handler_t)(int);
}

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
     // should this be 'volatile' ???
    static bool interruptFlag;
      //this is a pointer function that will point back to the original
      //signal handler if there was one.
    static msq_sig_handler_t oldSigIntHandler;
    
    
  };
  
  
}

#endif
