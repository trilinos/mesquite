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
/*!
  \file   MesquiteError.hpp
  \brief  Contains the error mechanism used in Mesquite.

  The MsqError class provides the error objects that are passed as the
  last argument of most Mesquite functions.
  The Error Mechanism includes:
  - stack tracing
  - context information (i.e. non standard error messages) 

  \author Thomas Leurent
  \date   2002-04-17
*/


#ifndef MesquiteError_hpp
#define MesquiteError_hpp

#include <new>
#include "Mesquite.hpp"
#include <string>
MSQ_USE(string);

/*! \def __SDIR__
  \brief Source file directory.
  
  Defines the directory where the source file is located; used
  in printing error messages. Each makefile has an entry
  LOCDIR     =  thedirectory
  and bmake/common includes in CFLAGS -D__SDIR__='"${LOCDIR}"'
  which is a flag passed to the compilers.
*/
#if !defined(__SDIR__)
#define __SDIR__ "unknowndirectory"
#endif
 
/*! \def __FUNC__
  \brief Redefined with the function name before each function.
  
  This is part of the stack tracing mechanism. __FUNC__ must be redefined
  with the function name before each function that
  uses the MSQ_CHKERR(err) macro. __FUNC__ is used
  in printing error messages.
*/
#if !defined(__FUNC__)
#define __FUNC__ "unknownfunction"
#endif

namespace Mesquite
{

  void out_of_store();
  
/*! \enum Error_Codes
  \brief Generic error codes.
  
  These error codes are used in
  many different places in the MESQUITE source code.
*/
enum Error_Codes {
  MSQ_NO_ERROR           =  0,
  MSQ_PRINT_STACK        = -1,
  MSQ_DO_NOT_PRINT       = -2,
  MSQ_MEM_ERR            = 55,   /* unable to allocate the requested memory */
  MSQ_NULL_ERR           = 56,   /* null data pointer */
  MSQ_INPUT_ERR          = 57,   /* something is wrong with input to function */
  MSQ_INIT_ERR           = 58,   /* data structure not initialized */
  MSQ_FILE_OPEN_ERR      = 59,   /* unable to open file */
  MSQ_FREE_ERR           = 60,   /* unable to free memory */
  MSQ_INVALID_MESH_ERR   = 61,   /* unable to free memory */
  MSQ_DIVIDE_BY_ZERO_ERR = 62,   /* division by zero */
  MSQ_DATA_ERR           = 63,    /* incorrect data */
  MSQ_NO_PD_STORAGE_MODE = 100    /* no storage mode chosen within PatchData */
};


/*! \class MsqError
  \brief Mesquite's error object.

  A MsqError object is passed by reference to all the functions that
  use the error handling mechanism. In order to obtain succesful stack
  tracing, a function must have an MsqError argument if it calls a
  function which has an MsqError argument itself.

  An MsqError object can be activated in 3 ways, which can be combined:
  -# err.errorOn = true;
  -# err.set_msg("bad idea");
  -# err.set_error_code(MSQ_MEM_ERR);
  */
class MsqError {
public:
  // Constructors
  MsqError() : errorOn(false), printError(true), errorCode(MSQ_NO_ERROR), msg("") {} // No error by default
  //! Maybe this constructor isn't useful.

  MsqError(bool , enum Error_Codes, string = "");

  //! Maybe this constructor isn't useful.
  MsqError(Error_Codes, string = "");
  
  //! sets a context dependent (non-generic) error message.
  void set_msg(const string&);  // also sets errorOn to true
  //! sets a generic error code, later handled by MsqError::handler.
  void set_error_code(enum Error_Codes); // also sets errorOn to true

  //! resets error object to non-active state (no error).
  void reset(); 
  
  //! Error Handler - Accessed through MSQ_CHKERR(err)
  void  handler(int line, const char *func,
                const char* file, const char *dir);
  
  // data member with direct access
  bool errorOn;

    /*! If printError is true, then the next time the error object is
      checked and an error is on, we print an error message.  If printError
      is false, we do not print the message.
      If MESQUITE_PRINT_ERROR_STACK is defined then printError
      is always true.  Otherwise, printError is only true if no error
      messages have been printed from this object (at least since the
      last time the object was reset).
    */
  bool printError;

protected:
    //! Copy constructor -- protected access
  MsqError(const MsqError &A);

private:
  enum Error_Codes errorCode;  
  string msg;
};

} // namespace


/*! \def MSQ_CHKERR(err)
  \brief Mesquite's Error Checking macro.

  This macro should be used after a call to any function
  that contains a MsqError (trailing) argument.
  Placing the macro on the same lien as the function will
  allow to print at runtime the exact line on which the error occured.  
*/
#define MSQ_CHKERR(err)     {if (err.errorOn) err.handler(__LINE__,__FUNC__,__FILE__,__SDIR__);}

/*! \def MSQ_CHECK_NULL(a)
  \brief checks for null integer, pointer, etc...
*/
#ifndef MSQ_CHECK_NULL
#define MSQ_CHECK_NULL(a) \
{ \
   if (a==0) { MsqError e(true, MSQ_NULL_ERR); MSQ_CHKERR(e); }\
}
#endif

#endif


