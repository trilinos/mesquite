#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
#include <vector>

//delete the following includes whenever we move
//print_timing_diagnostics()
#include <iostream>
#include <iomanip>
// Forward declarations
static void default_print_info_func(const char* msg);
static void default_print_warning_func(const char* msg);
static void default_print_error_func(const char* msg);

typedef void (*MsqInfoFunc)(const char* msg);

// Static member variables - the function pointers
MsqInfoFunc Mesquite::Message::printInfoFunc =
  default_print_info_func;
MsqInfoFunc Mesquite::Message::printWarningFunc =
  default_print_warning_func;
MsqInfoFunc Mesquite::Message::printErrorFunc =
  default_print_error_func;

// Functions to set the printing function pointers
void Mesquite::Message::set_print_info(void (*print_func)(const char*))
{ printInfoFunc = print_func; }

void Mesquite::Message::set_print_warning(void (*print_func)(const char*))
{ printWarningFunc = print_func; }

void Mesquite::Message::set_print_error(void (*print_func)(const char*))
{ printErrorFunc = print_func; }

// Default printing functions
static void default_print_info_func(const char* msg)
{ fputs(msg, stdout); }

static void default_print_warning_func(const char* msg)
{ fputs(msg, stdout); }

static void default_print_error_func(const char* msg)
{ fputs(msg, stderr); }

/*!Print a 'report' of the Timing data for StopWatches in the global
  StopWatchCollection.  */
void Mesquite::Message::print_timing_diagnostics()
{
  std::vector<Mesquite::StopWatchCollection::Key> sorted_keys;
  Mesquite::GlobalStopWatches.get_keys_sorted_by_time(sorted_keys);
  int number_of_keys=sorted_keys.size();
  int i =0;
  std::cout<<"\nTIME     | NUM. STARTS | TIMER NAME ("<<number_of_keys<<" timers)\n";
  for(i=0;i<number_of_keys;++i){
    std::cout<<std::setiosflags(std::ios::left)
             <<std::setw(10)
             <<Mesquite::GlobalStopWatches.total_time(sorted_keys[i])
             <<" "
             <<std::setw(13)
             <<Mesquite::GlobalStopWatches.number_of_starts(sorted_keys[i])
             <<" "
             <<Mesquite::GlobalStopWatches.get_string(sorted_keys[i])
             <<std::endl;
  }
}


