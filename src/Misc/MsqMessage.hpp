//- Class: MsqMessage
//- Description: MsqMessage class - used for reporting messages to
//-              the user.
//- Owner: Darryl Melander
//- Checked By:
//- Version:

#ifndef MSQMESSAGE_HPP
#define MSQMESSAGE_HPP

#ifdef USE_C_PREFIX_INCLUDES
#include <cstdio>
#include <cstdarg>
#else
#include <stdio.h>
#include <stdarg.h>
#endif

#define PRINT_INFO Mesquite::Message::print_info
#define PRINT_WARNING Mesquite::Message::print_warning
#define PRINT_ERROR Mesquite::Message::print_error

namespace Mesquite
{
  class Message
  {
  public:
    static void print_info(const char* format, ...);
    static void print_warning(const char* format, ...);
    static void print_error(const char* format, ...);
    
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

