// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 bu Tom Leurent
//  LAST-MOD: 13-Jan-03 at 11:32:08 by Todd Munson
//
// DESCRIPTION:
// ============
/*! \file MsqHessian.hpp

describe MsqHessian.hpp here

 */


#ifndef MsqHessian_hpp
#define MsqHessian_hpp

#include "Mesquite.hpp"
#include "Matrix3D.hpp"
#include "PatchData.hpp"

#ifdef USE_C_PREFIX_INCLUDES
#include <cassert>
#else
#include <assert.h>
#endif

#include <iostream>


namespace Mesquite
{
  
  /*!
    \class MsqHessian
    \brief Vector3D is the object that effeciently stores the objective function
    Hessian each entry is a Matrix3D object (i.e. a vertex Hessian). 
  */
  class MsqHessian
  {
  private:
    Matrix3D* mEntries;	   // size: number of nonzero blocks 
    size_t* mColIndex;     //< column indexes of the entries in the row. 
    size_t* mColInstr;	   //< accumulation pattern instructions
    size_t* mRowStart;	   // size: number of vertices

    
  public:
    void initialize(PatchData &pd, MsqError &err);
  };

} // namespace

#endif // MsqHessian_hpp
