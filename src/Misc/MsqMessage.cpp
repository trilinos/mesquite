#include "MsqMessage.hpp"

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
