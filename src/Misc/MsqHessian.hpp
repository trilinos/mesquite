// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 bu Thomas Leurent
//  LAST-MOD: 23-Jan-03 at 16:21:09 by Thomas Leurent
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
    PatchData* origin_pd;
    
    Matrix3D* mEntries;	   //!< CSR block entries.  size: number of nonzero blocks 
    size_t* mRowStart;	   //!< start of each row in mEntries. size: number of vertices.
    size_t* mColIndex;     //!< CSR block structure: column indexes of the row entries. 

    int* mAccumulation;	   //!< accumulation pattern instructions
    size_t* mAccumElemStart;  //!< Starting index in mAccumulation for element i, i=1,...

    int mSize; //!< number of rows (or number of columns, this is a square matrix).

  public:
    MsqHessian();
    
    void initialize(PatchData &pd, MsqError &err);
    int size() {return mSize;}
    //! returns the diagonal blocks, memory must be allocated before call.
    void get_diagonal_blocks(std::vector<Matrix3D> &diag, MsqError &err);
    void accumulate_entries(PatchData &pd, size_t elem_index,
                            Matrix3D mat3d_array[], size_t nb_mat3d, MsqError &err); 
    //! Hessian - vector product summed with a second vector.
    friend void axpy(Vector3D res[], int size_r,
                     const MsqHessian &H, const Vector3D x[], int size_x,
                     const Vector3D y[], int size_y, MsqError &err);
  };

  
#undef __FUNC__
#define __FUNC__ "axpy"
  /*!
    \param res: array of Vector3D in which the result is stored.
    \param x: vector multiplied by the Hessian.
    \param y: vector added to the Hessian vector product.
  */
  inline void axpy(Vector3D res[], int size_r,
                   const MsqHessian &H, const Vector3D x[], int size_x,
                   const Vector3D y[], int size_y, MsqError &err)
  {
    if ((size_r != H.mSize) || (size_x != H.mSize) || (size_y != H.mSize)) {
      // throw an error
    }

    Vector3D xl, yl, ml;

    const int nn = H.mSize;
    int rl; // row length
    int e;  // entries index
    int c;  // column index
    int i, j;
    size_t* col = H.mColIndex;
  
    for(i=0; i<nn; ++i)
      res[i] = y[i];
   
    for (i = 0; i < nn; ++i) {
      rl = H.mRowStart[i+1] - H.mRowStart[i];
      e = H.mRowStart[i];
      if (rl) {
       
        // Diagonal entry
        res[i] += H.mEntries[e]*x[i]; 
        ++e;
        assert(*col == i);

        //Non-diagonal entries
        for (j = 1; j < rl; ++j) {
          c = *(++col);
          res[i] += H.mEntries[e] * x[c];
          res[c] += transpose(H.mEntries[e]) * x[i];
          ++e;
        }
        ++col;
      }
    }
  }

} // namespace

#endif // MsqHessian_hpp
