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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

#include "MesquiteInterrupt.hpp"

bool Mesquite::MesquiteInterrupt::interruptFlag = false;
msq_sig_handler_t Mesquite::MesquiteInterrupt::oldSigIntHandler = SIG_DFL;


#ifdef ENABLE_INTERRUPT

//function that catches the interrupt
extern "C" void msq_sigint_handler(int)
{
    //notify the user that we recognized the interrupt signal
  Mesquite::Message::print_info("Handling interrupt signal.\n");
    //if there was an original handler, call it here.
  if( Mesquite::MesquiteInterrupt::oldSigIntHandler != SIG_DFL &&
      Mesquite::MesquiteInterrupt::oldSigIntHandler != SIG_ERR ){
    Mesquite::Message::print_info("Calling original signal handler\n");
    Mesquite::MesquiteInterrupt::oldSigIntHandler(SIGINT);
  }
    //set the signal handler again
  if( signal( SIGINT, msq_sigint_handler ) == SIG_ERR )
    Mesquite::Message::print_info("Can't catch SIGINT!\n");
    //set the flag that will be checked elsewhere in the code
  Mesquite::MesquiteInterrupt::interruptFlag = true;
  
}

// Set signal handler for SIGINT.  Either catch the interrupt or don't.
void Mesquite::MesquiteInterrupt::catch_interrupt( bool yesno )
{
  if( yesno )
  {
      //also save the pointer to the original handler
    Mesquite::MesquiteInterrupt::oldSigIntHandler = signal( SIGINT, msq_sigint_handler );
    if( Mesquite::MesquiteInterrupt::oldSigIntHandler == SIG_ERR )
      Mesquite::Message::print_info("Can't catch SIGINT!\n");
  }
  else
  {
      //make sure we are not sending SIG_ERR to signal
    if(Mesquite::MesquiteInterrupt::oldSigIntHandler == SIG_ERR )
      Mesquite::MesquiteInterrupt::oldSigIntHandler = SIG_DFL;
      //reset the signal handler.
    if( signal( SIGINT, Mesquite::MesquiteInterrupt::oldSigIntHandler )
        == SIG_ERR )
      Mesquite::Message::print_info("Can't reset SIGINT handler!\n");
  }
}
#endif
