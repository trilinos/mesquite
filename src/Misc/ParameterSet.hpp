#ifndef MESQUITE_PARAMETER_SET_HPP
#define MESQUITE_PARAMETER_SET_HPP


#ifdef USE_STD_INCLUDES
#include <cstddef>
#else
#include <stddef.h>
#endif

#include "Mesquite.hpp"

namespace Mesquite
{
  
  class ParameterSet
  {
  public:
    ParameterSet();
    ~ParameterSet();
    
    void add_int_parameter(const char* name,
                           int initial_value, MsqError &err);
    void set_int_parameter(const char* name,
                           int value, MsqError &err);
    void get_int_parameter(const char* name,
                           int* value, MsqError &err);
    
    void remove_parameter(const char* name, MsqError &err);
    
  private:
    union ParameterValue
    {
      int intVal;
      double dblVal;
      void* ptrVal;
      char* strVal;
      bool boolVal;
    };

    enum ParameterType
    {
      MSQ_INT,
      MSQ_DBL,
      MSQ_PTR,
      MSQ_STRING,
      MSQ_BOOL
    };
    
    struct ParameterRecord
    {
      char* name;
      ParameterType type;
      ParameterValue value;
    };
    
    ParameterRecord* mParameterArray;
    size_t mNumParameters;
    
      // returns the 0-based index of where the parameter
      // with the given name can be found in mParameterArray,
      // or mNumParameters if it can't be found.
    size_t get_parameter_index(const char* name);
    
    void generic_add_parameter(const char* name, MsqError &err);
  };

}

#endif
