// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

#include "MesquiteInterrupt.hpp"

bool Mesquite::MesquiteInterrupt::interruptFlag = false;

#ifdef ENABLE_INTERRUPT

//function that catches the interrupt
extern "C" void msq_sigint_handler(int)
{
  if( signal( SIGINT, msq_sigint_handler ) == SIG_ERR )
    PRINT_INFO("Can't catch SIGINT!\n");
  Mesquite::MesquiteInterrupt::interruptFlag = true;
  PRINT_INFO("Handling interrupt signal.\n");
}

// Set signal handler for SIGINT.  Either catch the interrupt or don't.
void Mesquite::MesquiteInterrupt::catch_interrupt( bool yesno )
{
  if( yesno )
  {
    if( signal( SIGINT, msq_sigint_handler ) == SIG_ERR )
      PRINT_INFO("Can't catch SIGINT!\n");
  }
  else
  {
    if( signal( SIGINT, SIG_DFL ) == SIG_ERR )
      PRINT_INFO("Can't reset SIGINT handler!\n");
  }
}
#endif
