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
#ifndef MESQUITE_HPP
#define MESQUITE_HPP

#ifdef WIN32
#pragma warning ( 4 : 4786)
#endif

#include <stdexcept>
#ifdef USE_STD_INCLUDES
#include <iostream>
#else
#include <iostream.h>
#endif

#ifdef USE_C_PREFIX_INCLUDES
#include <cmath>
#include <cstddef>
#include <cstring>
#else
#include <math.h>
#include <stddef.h>
#include <string.h>
#endif

#include <float.h>

/*! \file Mesquite.hpp
 */

/*!
  \namespace Mesquite
  Copyright 2003 Sandia Corporation and the University of Chicago. Under
  the terms of Contract DE-AC04-94AL85000 with Sandia Corporation and
  Contract W-31-109-ENG-38 with the University of Chicago, the U.S.
  Government retains certain rights in this software.

*/
namespace Mesquite
{
//#ifndef MESQUITE_PRINT_ERROR_STACK
//#define MESQUITE_PRINT_ERROR_STACK
//#endif
#ifndef ENABLE_INTERRUPT
#define ENABLE_INTERRUPT
#endif  
  typedef int StatusCode;

  typedef double real;

  enum StatusCodeValues
  {
    MSQ_FAILURE = 0,
    MSQ_SUCCESS
  };

  enum EntityTopology
  {
    POLYGON =7,
    TRIANGLE =8,
    QUADRILATERAL =9,
    POLYHEDRON =10,
    TETRAHEDRON =11,
    HEXAHEDRON =12,
    PRISM =13,
    PYRAMID =14,
    SEPTAHEDRON =15,
    MIXED
  };

    // Version information
  const char* version_string(bool include_build_number = false);
  unsigned int major_version_number();
  unsigned int minor_version_number();
  unsigned int build_number();
  enum ReleaseType
  {
    STABLE_RELEASE,
    BETA,
    ALPHA
  };
  Mesquite::ReleaseType release_type();
  
    // This function should be called by the calling
    // application when it wants to interrupt a Mesquite
    // algorithm before it has completed, such as when
    // the user hits ctrl-c
  void signal_interrupt();
  
    //GLOBAL variables
  const int MSQ_MAX_NUM_VERT_PER_ENT=8;
  const int MSQ_HIST_SIZE=7;//number of division in histogram
  static const double MSQ_SQRT_TWO = sqrt(2.0);
  static const double MSQ_SQRT_THREE = sqrt(3.0);
  static const double MSQ_SQRT_THREE_DIV_TWO=MSQ_SQRT_THREE/2.0;
  static const double MSQ_SQRT_THREE_INV=1.0/MSQ_SQRT_THREE;
  static const double MSQ_SQRT_TWO_INV=1.0/MSQ_SQRT_TWO;
  static const double MSQ_SQRT_TWO_DIV_SQRT_THREE=MSQ_SQRT_TWO/MSQ_SQRT_THREE;
  static const double MSQ_ONE_THIRD = 1.0 / 3.0;
  static const double MSQ_TWO_THIRDS = 2.0 / 3.0;

#ifdef INT_MAX
  const int MSQ_INT_MAX = INT_MAX;
#else     
  const int MSQ_INT_MAX = 2147483647;
#endif

#ifdef INT_MIN
const int MSQ_INT_MIN = INT_MIN;
#else
const int MSQ_INT_MIN = -2147483647;
#endif

#ifdef DBL_MIN
  const double MSQ_DBL_MIN = DBL_MIN;
#else
  const double MSQ_DBL_MIN = 1.0E-30;
#endif
  const double MSQ_MIN=DBL_MIN;
  
#ifdef DBL_MAX
  const double MSQ_DBL_MAX = DBL_MAX;
#else
  const double MSQ_DBL_MAX = 1.0E30;
#endif    
  const double MSQ_MAX = DBL_MAX;
  const double MSQ_MAX_CAP = 1.e6;

    //macro to return the min/max of a set of arguements.  The integer
    // (e.g., '2') tells how many arguements should be passed for comparison.
#ifndef MSQ_MIN_2 
#define MSQ_MIN_2(a,b)     ( (a) < (b) ? (a) : (b) )
#endif
#ifndef MSQ_MAX_2           
#define MSQ_MAX_2(a,b)     ( (a) > (b) ? (a) : (b) )
#endif
}

#ifndef USE_FUNCTION_TIMERS
#define USE_FUNCTION_TIMERS
#endif

#ifdef USE_STD_INCLUDES
  #define MSQ_USE(A) using std::A 
#else
  #define MSQ_USE(A) //empty statement
#endif


#define MSQ_DEBUG
#define MSQ_DEBUG4 // notify of all constructor / copy / destructor uses


/* the Mesquite Debugging system
     Level 0 provides no information and the debug macros are empty
     Level 1 provides user function information only
             e.g. threshold set, function used, etc.
     Level 2 provides basic algorithmic information for each local
             submesh
     Level 3 provides more information, data structures, and details
             than most users would want to know about 

     The default is Level 0
*/
//#define MSQ_DBG3

#ifdef MSQ_DBG3
#define MSQ_DEBUG_LEVEL 3
#elif  MSQ_DBG2
#define MSQ_DEBUG_LEVEL 2
#elif  MSQ_DBG1
#define MSQ_DEBUG_LEVEL 1
#elif  MSQ_DBG0
#define MSQ_DEBUG_LEVEL 0
#else
#define MSQ_DEBUG_LEVEL 0
#endif

#if MSQ_DEBUG_LEVEL != 0
#define MSQ_DEBUG_PRINT(level, statement)\
{\
   if ((level <= MSQ_DEBUG_LEVEL))\
     {\
     fprintf(stdout,statement);\
     fflush(stdout);\
     }\
}
#define MSQ_DEBUG_ACTION(level, action)\
{\
   if ((level <= MSQ_DEBUG_LEVEL))\
     {\
     action\
     fflush(stdout);\
     }\
}
#else
#define MSQ_DEBUG_PRINT(level, statement)\
{\
}
#define MSQ_DEBUG_ACTION(level, statement)\
{\
}
#endif

#endif
