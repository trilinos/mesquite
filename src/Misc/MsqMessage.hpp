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
//- Class: MsqMessage
//- Description: MsqMessage class - used for reporting messages to
//-              the user.
//- Owner: Darryl Melander
//- Checked By:
//- Version:

#ifndef MSQMESSAGE_HPP
#define MSQMESSAGE_HPP

/* This code is broken.  The difference between the no-extension and .h
   headers is that the former define all standard C functions in the
   "std" namespace.  The no-extension headers are included but the
   namespace is not specified when the corresponding symbols are used.
   
   Also, it is more portable to just use the .h headers all the time.
   Do that, at least for now. -- J.Kraftcheck
   FIXME
#ifdef USE_C_PREFIX_INCLUDES
#include <cstdio>
#include <cstdarg>
#else
*/
#include <stdio.h>
#include <stdarg.h>
/*
#endif
*/


namespace Mesquite
{
  class Message
  {
  public:
    static void print_info(const char* format, ...);
    static void print_warning(const char* format, ...);
    static void print_error(const char* format, ...);

    static void print_timing_diagnostics();
    
    static void set_print_info(void (*print_func)(const char*));
    static void set_print_warning(void (*print_func)(const char*));
    static void set_print_error(void (*print_func)(const char*));
    
  private:
    static void (*printInfoFunc)(const char* msg);
    static void (*printWarningFunc)(const char* msg);
    static void (*printErrorFunc)(const char* msg);
  };
  
  inline void Message::print_info(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    
    char msgbuf[10240];
    vsprintf(msgbuf, format, args);
    
    (*printInfoFunc)(msgbuf);
    
    va_end(args);
  }
  
  inline void Message::print_error(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    
    char msgbuf[10240];
    vsprintf(msgbuf, format, args);
    
    (*printErrorFunc)(msgbuf);
    
    va_end(args);
  }
  
  inline void Message::print_warning(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    
    char msgbuf[10240];
    vsprintf(msgbuf, format, args);
    
    (*printWarningFunc)(msgbuf);
    
    va_end(args);
  }
}

#endif

