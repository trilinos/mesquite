// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 bu Thomas Leurent
//  LAST-MOD: 21-Mar-03 at 11:19:14 by Thomas Leurent
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
  class ObjectiveFunction;
  
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
    
    Matrix3D* mPreconditionner;
    int precondArraySize;
    
    Vector3D* r; //!< array used in the CG solver
    Vector3D* z; //!< array used in the CG solver
    Vector3D* p; //!< array used in the CG solver
    Vector3D* w; //!< array used in the CG solver
    int cgArraySizes; //!< size of arrays allocated in the CG solver.
    int maxCGiter; //!< max nb of iterations of the CG solver.
    
  public:
    MsqHessian();
    ~MsqHessian();
    
    void initialize(PatchData &pd, MsqError &err);
    int size() {return mSize;}
    //! returns the diagonal blocks, memory must be allocated before call.
    void get_diagonal_blocks(std::vector<Matrix3D> &diag, MsqError &err);
    Matrix3D* get_block(size_t i, size_t j);
    void accumulate_entries(PatchData &pd, size_t elem_index,
                            Matrix3D mat3d_array[], MsqError &err);
    void compute_preconditionner(MsqError &err);
    void apply_preconditionner(Vector3D zloc[], Vector3D rloc[], MsqError &err);
    void cg_solver(Vector3D x[], Vector3D b[], MsqError &err);
    //! Hessian - vector product, summed with a second vector (optional).
    friend void axpy(Vector3D res[], int size_r,
                     const MsqHessian &H, const Vector3D x[], int size_x,
                     const Vector3D y[], int size_y, MsqError &err);
    friend class ObjectiveFunction;
  };

  
#undef __FUNC__
#define __FUNC__ "axpy"
  /*!
    \param res: array of Vector3D in which the result is stored.
    \param size_r: size of the res array.
    \param x: vector multiplied by the Hessian.
    \param size_x: size of the x array.
    \param y: vector added to the Hessian vector product. Set to 0 if not needed.
    \param size_y: size of the y array. Set to 0 if not needed.
  */
  inline void axpy(Vector3D res[], int size_r,
                   const MsqHessian &H, const Vector3D x[], int size_x,
                   const Vector3D y[], int size_y, MsqError &err)
  {
    if ((size_r != H.mSize) || (size_x != H.mSize) ||
        (size_y != H.mSize && size_y != 0)) {
      // throw an error
    }

    const int nn = H.mSize;
    int rl; // row length
    int e;  // entries index
    int c;  // column index
    int i, j;
    size_t* col = H.mColIndex;

    if (y!=0) 
      for(i=0; i<nn; ++i)
        res[i] = y[i];
    else // y==0
      for(i=0; i<nn; ++i)
        res[i] = 0.;
   
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



#undef __FUNC__
#define __FUNC__ "MsqHessian::apply_preconditionner"
  /*! Computes \f$ z=M^{-1}r \f$ . */
  inline void MsqHessian::apply_preconditionner(Vector3D zloc[],
                                                Vector3D rloc[],
                                                MsqError& /*err*/)
  {
    int m;
    // preconditionner is identity matrix for now.
    for (m=0; m<mSize; ++m) {
      zloc[m] = mPreconditionner[m] * rloc[m]; 
    }
  }

  
#undef __FUNC__
#define __FUNC__ "MsqHessian::get_block"
  /*! Computes \f$ z=M^{-1}r \f$ . */
  inline Matrix3D* MsqHessian::get_block(size_t i, size_t j)
  {
    assert(i<mSize);
    assert(j<mSize);
    assert(j>=i);

//    std::cout << "mColIndex[mRowStart[i]+j]: " << mColIndex[mRowStart[i]+j] << "  mRowStart[i]: " <<mRowStart[i]<< " for i=" << i << " j=" << j << std::endl; //dbg
    assert( mColIndex[mRowStart[i]+j] == j ); // since our Matrix is in fact not sparse. 
    
    return ( mEntries + (mRowStart[i]+j)  );
  }

  
} // namespace

#endif // MsqHessian_hpp
