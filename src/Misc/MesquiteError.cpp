#include "MesquiteError.hpp"

using namespace Mesquite;

MSQ_USE(cerr);
MSQ_USE(endl);

/* Special handler for the bad_alloc exception thrown by new */ 
void Mesquite::out_of_store()
{
  cerr << "MESQUITE ERROR: Operator new failed: out of store.\n";
#if USE_STD_INCLUDES
  throw std::bad_alloc();
#else
  throw bad_alloc();
#endif  
}

// Constructors
MsqError::MsqError(bool e, enum Error_Codes ec, string m)
{
  errorOn = e;
  printError=true;
  errorCode = ec;
  if (m!="")
    msg = m;
}

MsqError::MsqError(Error_Codes ec, string m)
{
  errorOn = true;
  printError=true;
  errorCode = ec;
  if (m != string(""))
    msg = m;
}

/*! 
The copy constructor has protected access since users should always pass
errors by reference.
The copy constructor sends a warning if used, since MsqError object
should always be passed by reference.
*/
MsqError::MsqError(const MsqError &A)
{
  cerr << "WARNING: the MsqError() copy constructor has been used\n"
            << "         MsqError objects should always be passed by reference";
  errorOn = A.errorOn;
  errorCode = A.errorCode;
  if (A.msg!="")
    msg = A.msg;
  printError=A.printError;
}


void MsqError::set_msg(const string &m)
{
  errorOn = true;
  msg = m;
}

void MsqError::set_error_code(enum Error_Codes ec)
{
  errorOn = true;
  errorCode = ec;
}

void MsqError::reset()
{
  errorOn = false;
  errorCode = MSQ_NO_ERROR;
  msg = "";
  printError=true;
}

/*! \fn MsqError::handler(int line, const char *func,const char* file, const char *dir)
  This functions provides by default the following error handling mechanism:
  - prints the contextual error message, if any was set with set_msg(...) .
  - handles a standard error code, if any was set with set_error_code(...) .
  - prints stack tracing information.
  */
#undef __FUNC__
#define __FUNC__ "MsqError::handler"
void MsqError::handler(int line, const char *func,
                       const char* file, const char *dir)
{
  
  if(printError){
    if (msg!="" && errorCode!=MSQ_PRINT_STACK) {
      cerr << "MESQUITE ERROR: " << msg << endl;
    }
    switch(errorCode){
      case MSQ_MEM_ERR:
        cerr << "MESQUITE ERROR:  Out of memory. \n";
        break;
      case MSQ_NULL_ERR:
        cerr << "MESQUITE ERROR:  Null pointer. \n";
        break;
      case MSQ_INIT_ERR:
        cerr << "MESQUITE ERROR:  Data Structure Not Initialized. \n";
        break;
      case MSQ_INPUT_ERR:
        cerr << "MESQUITE ERROR:  Incorrect Input \n";
        break;
      case MSQ_FILE_OPEN_ERR:
        cerr << "MESQUITE ERROR:  File open error \n";
        break;
      case MSQ_FREE_ERR:
        cerr << "MESQUITE ERROR:  Error freeing memory \n";
        break;
      case MSQ_INVALID_MESH_ERR:
        cerr << "MESQUITE ERROR:  Invalid Mesh; use SMuntangle to create a valid mesh prior to smoothing \n";
        break;
      case MSQ_DIVIDE_BY_ZERO_ERR:
        cerr << "MESQUITE ERROR:  Division by zero \n";
        break;
      case MSQ_DATA_ERR:
        cerr << "MESQUITE ERROR:  Incorrect data \n";
        break;
      case MSQ_NO_PD_STORAGE_MODE:
        cerr << "MESQUITE ERROR: no storage mode set in PatchData object.\n";
        break;
      case MSQ_NO_ERROR:
        break;
      default:
        break;
        
    }
    cerr << "MESQUITE ERROR: " << func << "()  line " << line
              << " in " << dir <<"/"<< file << endl;
  
#ifdef MESQUITE_PRINT_ERROR_STACK
    errorCode = MSQ_PRINT_STACK; /* set it to print the stack */
#else
    errorCode = MSQ_DO_NOT_PRINT;
    printError=false;
#endif
  }
  
  return ;
}

