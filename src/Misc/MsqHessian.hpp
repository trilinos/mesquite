// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 bu Thomas Leurent
//  LAST-MOD: 17-Jan-03 at 17:25:29 by Thomas Leurent
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
  protected:
    Matrix3D* mEntries;	   //!< CSR block entries.  size: number of nonzero blocks 
    size_t* mRowStart;	   //!< start of each row in mEntries. size: number of vertices.
    size_t* mColIndex;     //!< CSR block structure: column indexes of the row entries. 

    int* mAccumulation;	   //!< accumulation pattern instructions

    int mSize; //!< number of rows (or number of columns, this is a square matrix).
    
  public:
    void initialize(PatchData &pd, MsqError &err);
    int size() {return mSize;}
    //! returns the diagonal blocks, memory must be allocated before call.
    void get_diagonal_blocks(std::vector<Matrix3D> &diag, MsqError &err);
  };

} // namespace

#endif // MsqHessian_hpp
