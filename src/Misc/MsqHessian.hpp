// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 bu Thomas Leurent
//  LAST-MOD:  9-Apr-03 at 17:11:27 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqHessian.hpp

 The MsqHessian class stores a sparse hessian for a given objective 
 function. The objective function must be C2 and such that its hessian
 has non-zero entries only for the duplet of derivatives corresponding 
 to nodes of a same element. 

 \author Thomas Leurent
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
    
    Matrix3D* mEntries;	//!< CSR block entries. size: nb of nonzero blocks, i.e. mRowStart[mSize] . 
    size_t* mRowStart;	//!< start of each row in mEntries. size: nb of vertices (mSize).
    size_t* mColIndex;  //!< CSR block structure: column indexes of the row entries. 

    int* mAccumulation;	   //!< accumulation pattern instructions
    size_t* mAccumElemStart;  //!< Starting index in mAccumulation for element i, i=1,...

    size_t mSize; //!< number of rows (or number of columns, this is a square matrix).
    
    Matrix3D* mPreconditioner;
    size_t precondArraySize;
    
    Vector3D* r; //!< array used in the CG solver
    Vector3D* z; //!< array used in the CG solver
    Vector3D* p; //!< array used in the CG solver
    Vector3D* w; //!< array used in the CG solver
    size_t cgArraySizes; //!< size of arrays allocated in the CG solver.
    size_t maxCGiter; //!< max nb of iterations of the CG solver.
    
  public:
    MsqHessian();
    ~MsqHessian();
    
    void initialize(PatchData &pd, MsqError &err);
    void zero_out();
    size_t size() {return mSize;}
    //! returns the diagonal blocks, memory must be allocated before call.
    void get_diagonal_blocks(std::vector<Matrix3D> &diag, MsqError &err);
    Matrix3D* get_block(size_t i, size_t j);
    void accumulate_entries(PatchData &pd, size_t elem_index,
                            Matrix3D mat3d_array[], MsqError &err);
    void compute_preconditioner(MsqError &err);
    void apply_preconditioner(Vector3D zloc[], Vector3D rloc[], MsqError &err);
    void cg_solver(Vector3D x[], Vector3D b[], MsqError &err);
    //! Hessian - vector product, summed with a second vector (optional).
    friend void axpy(Vector3D res[], size_t size_r,
                     const MsqHessian &H, const Vector3D x[], size_t size_x,
                     const Vector3D y[], size_t size_y, MsqError &err);
    friend class ObjectiveFunction;
    friend std::ostream& operator<<(std::ostream &s, const MsqHessian &h);
  };


#undef __FUNC__
#define __FUNC__ "MsqHessian::zero_out"
  /*! Sets all Hessian entries to zero. This is usually used before 
      starting to accumulate elements hessian in the objective function
      hessian. */
  inline void MsqHessian::zero_out()
  {
    if (mSize==0) return; // empty hessian.
    
    size_t i;
    for (i=0; i<mRowStart[mSize]; ++i) {
      mEntries[i] = 0.;
    }
  }

  
#undef __FUNC__
#define __FUNC__ "MsqHessian::accumulate_entries"
  /*! Accumulates entries of an element hessian into an objective function
      hessian. Make sure to use zero_out() before starting the accumulation
      process. 

    \param pd: PatchData in that contains the element which Hessian
    we are accumulating in the Hessian matrix. This must be the same
    PatchData that was used in MsqHessian::initialize().
    \param elem_index: index of the element in the PatchData.
    \param mat3d_array: This is the upper triangular part of the element Hessian 
    for all nodes, including fixed nodes, for which the entries must be null Matrix3Ds.
    \param nb_mat3d. The size of the mat3d_array: (n+1)n/2, where n is
    the number of nodes in the element.
*/
  inline void MsqHessian::accumulate_entries(PatchData &pd, size_t elem_index,
                                    Matrix3D mat3d_array[], MsqError &err)
  {
    MSQ_DEBUG_ACTION(1,{if (&pd != origin_pd) {
      err.set_msg("Cannot accumulate elements from a different patch. "
                  "Use MsqHessian::initialize first."); return;}});

    size_t nve = pd.get_element_array(err)[elem_index].vertex_count(); MSQ_CHKERR(err);
    size_t nb_mat3d = (nve+1)*nve/2;
    
    size_t e = mAccumElemStart[elem_index];
    size_t i;
    for (i=0; i<nb_mat3d; ++i) {
       if (mAccumulation[e] >= 0)
          mEntries[mAccumulation[e]] += mat3d_array[i];
       else
          mEntries[-mAccumulation[e]].plus_transpose_equal(mat3d_array[i]);
       ++e;
    }
    
    assert( e == mAccumElemStart[elem_index+1] );
  }
   
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
  inline void axpy(Vector3D res[], size_t size_r,
                   const MsqHessian &H, const Vector3D x[], size_t size_x,
                   const Vector3D y[], size_t size_y, MsqError &err)
  {
    if ((size_r != H.mSize) || (size_x != H.mSize) ||
        (size_y != H.mSize && size_y != 0)) {
      // throw an error
    }

    const size_t nn = H.mSize;
    size_t rl; // row length
    size_t e;  // entries index
    size_t c;  // column index
    size_t i, j;
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
#define __FUNC__ "MsqHessian::apply_preconditioner"
  /*! Computes \f$ z=M^{-1}r \f$ . */
  inline void MsqHessian::apply_preconditioner(Vector3D zloc[],
                                                Vector3D rloc[],
                                                MsqError& /*err*/)
  {
    size_t m;

    for (m=0; m<mSize; ++m) {
#ifdef DIAGONAL_PRECONDITIONER
      // preconditioner is identity matrix for now.
      zloc[m][0] = mPreconditioner[m][0][0] * rloc[m][0]; 
      zloc[m][1] = mPreconditioner[m][1][1] * rloc[m][1]; 
      zloc[m][2] = mPreconditioner[m][2][2] * rloc[m][2]; 
#else
      // z = inv(L^T) * r
      zloc[m][0] = rloc[m][0];
      zloc[m][1] = rloc[m][1] - mPreconditioner[m][0][1] * zloc[m][0];
      zloc[m][2] = rloc[m][2] - mPreconditioner[m][0][2] * zloc[m][0] - mPreconditioner[m][1][2] * zloc[m][1];

      // z = inv(D) * z
      zloc[m][0] *= mPreconditioner[m][0][0];
      zloc[m][1] *= mPreconditioner[m][1][1];
      zloc[m][2] *= mPreconditioner[m][2][2];

      // z = inv(L) * z
      zloc[m][2] = zloc[m][2];
      zloc[m][1] = zloc[m][1] - mPreconditioner[m][1][2] * zloc[m][2];
      zloc[m][0] = zloc[m][0] - mPreconditioner[m][0][1] * zloc[m][1] - mPreconditioner[m][0][2] * zloc[m][2];
#endif
    }
  }

  
#undef __FUNC__
#define __FUNC__ "MsqHessian::get_block"
  /*! Returns a pointer to the Matrix3D block at position i,j if it exist. 
      Returns the NULL pointer if position i,j is a NULL entry.
      Note that block i,j must be in the upper triangular part of the 
      (symetric) hessian. */
  inline Matrix3D* MsqHessian::get_block(size_t i, size_t j)
  {
    size_t c;
    assert(i<mSize);
    assert(j<mSize);
    assert(j>=i);

    for (c=mRowStart[i]; c<mRowStart[i+1]; ++c) {
      if (mColIndex[c] == j)
        return ( mEntries + c );
    }
    
    // if there is no block at position i,j (zero entry).
    return NULL;
  }


/* ------------------ I/O ----------------- */

   //! Prints out the MsqHessian blocks.
  inline std::ostream& operator<<(std::ostream &s, const MsqHessian &h)
  {
    size_t i,j;
    s << "MsqHessian of size: " << h.mSize <<"x"<< h.mSize << "\n";
    for (i=0; i<h.mSize; ++i) {
      s << " ROW " << i << " ------------------------\n";
      for (j=h.mRowStart[i]; j<h.mRowStart[i+1]; ++j) {
        s << "   column " << h.mColIndex[j] << " ----\n";
        s << h.mEntries[j]; 
      } 
    }
    return s;
  }
  
} // namespace

#endif // MsqHessian_hpp
