// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 by Thomas Leurent
//  LAST-MOD:  9-Apr-03 at 09:27:54 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqHessian.cpp

describe MsqHessian.cpp here

*/


#include "MsqHessian.hpp"
#include "MsqTimer.hpp"

using namespace Mesquite;

using std::cout;
using std::endl;


MsqHessian::MsqHessian() :
  origin_pd(0), mEntries(0), mRowStart(0), mColIndex(0), 
  mAccumulation(0), mAccumElemStart(0), mSize(0), 
  mPreconditionner(0), precondArraySize(0),
  r(0), z(0), p(0), w(0), cgArraySizes(0), maxCGiter(50)
{ }


MsqHessian::~MsqHessian()
{
  delete[] mEntries;	
  delete[] mRowStart;	
  delete[] mColIndex; 

  delete[] mAccumulation;
  delete[] mAccumElemStart;

  delete[] mPreconditionner;

  delete[] r;
  delete[] z;
  delete[] p;
  delete[] w;
}

  
#undef __FUNC__
#define __FUNC__ "MsqHessian::initialize"
/*! \brief creates a sparse structure for a Hessian, based on the
  connectivity information contained in the PatchData.
  Only the upper triangular part of the Hessian is stored. */
void MsqHessian::initialize(PatchData &pd, MsqError &err)
{
  FUNCTION_TIMER_START(__FUNC__);
  delete[] mEntries;
  delete[] mRowStart;
  delete[] mColIndex;
  delete[] mAccumulation;
  delete[] mAccumElemStart;
  
  size_t num_vertices = pd.num_vertices();
  size_t num_elements = pd.num_elements();
  std::vector<size_t> vtx_list;
  size_t e, r, rs, re, c, cs, ce, nz, nnz, nve, i, j;
  MsqMeshEntity* element_array = pd.get_element_array(err); MSQ_CHKERR(err);

  if (num_vertices == 0) {
    err.set_msg("No vertices in PatchData");
    return;
  }

  mSize = num_vertices;

  // Calculate the offsets for a CSC representation of the accumulation
  // pattern.

  size_t* col_start = new size_t[num_vertices + 1];
  mAccumElemStart = new size_t[num_elements+1];
  mAccumElemStart[0] = 0;
  
  for (i = 0; i < num_vertices; ++i) {
    col_start[i] = 0;
  }

  for (e = 0; e < num_elements; ++e) {
    nve = element_array[e].vertex_count();
    element_array[e].get_vertex_indices(vtx_list);
    mAccumElemStart[e+1] = mAccumElemStart[e] + (nve+1)*nve/2;
    
    for (i = 0; i < nve; ++i) {
      r = vtx_list[i];
      
      for (j = i; j < nve; ++j) {
        c = vtx_list[j];

        if (r <= c) {
          col_start[c]++;
        }
        else {
          col_start[r]++;
        }
      }
    }
  }

  nz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = col_start[i];
    col_start[i] = nz;
    nz += j;
  }
  col_start[i] = nz;

  // Finished putting matrix into CSC representation

  int* row_instr = new int[5*nz];
  size_t* row_index = new size_t[nz];

  nz = 0;
  for (e = 0; e < num_elements; ++e) {
    nve = element_array[e].vertex_count();
    element_array[e].get_vertex_indices(vtx_list);

    for (i = 0; i < nve; ++i) {
      r = vtx_list[i];

      for (j = i; j < nve; ++j) {
        c = vtx_list[j];

        if (r <= c) {
          row_index[col_start[c]] = r;
          row_instr[col_start[c]] = nz;
          ++col_start[c];
        }
        else {
          row_index[col_start[r]] = c;
          row_instr[col_start[r]] = -nz;
          ++col_start[r];
        }
        
        ++nz;
      }
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    col_start[i+1] = col_start[i];
  }
  col_start[1] = col_start[0];
  col_start[0] = 0;

  //   cout << "col_start: ";
  //   for (int t=0; t<num_vertices+1; ++t)
  //     cout << col_start[t] << " ";
  //   cout << endl;
  //   cout << "row_index: ";
  //   for (int t=0; t<nz; ++t)
  //     cout << row_index[t] << " ";
  //   cout << endl;
  //   cout << "row_instr: ";
  //   for (int t=0; t<nz; ++t)
  //     cout << row_instr[t] << " ";
  //   cout << endl;
  
  
  // Convert CSC to CSR
  // First calculate the offsets in the row

  size_t* row_start = new size_t[num_vertices + 1];
    
  for (i = 0; i < num_vertices; ++i) {
    row_start[i] = 0;
  }

  for (i = 0; i < nz; ++i) {
    ++row_start[row_index[i]];
  }

  nz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = row_start[i];
    row_start[i] = nz;
    nz += j;
  }
  row_start[i] = nz;
    
  // Now calculate the pattern

  size_t* col_index = new size_t[nz];
  int* col_instr = new int[nz];

  for (i = 0; i < num_vertices; ++i) {
    cs = col_start[i];
    ce = col_start[i+1];

    while(cs < ce) {
      r = row_index[cs];

      col_index[row_start[r]] = i;
      col_instr[row_start[r]] = row_instr[cs];

      ++row_start[r];
      ++cs;
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    row_start[i+1] = row_start[i];
  }
  row_start[1] = row_start[0];
  row_start[0] = 0;

  delete[] row_index;

  // Now that the matrix is CSR
  // Column indices for each row are sorted

  // Compaction -- count the number of nonzeros
  mRowStart = col_start;   // don't need to reallocate
  mAccumulation = row_instr;   // don't need to reallocate

  for (i = 0; i <= num_vertices; ++i) {
    mRowStart[i] = 0;
  }

  nnz = 0;
  for (i = 0; i < num_vertices; ++i) {
    rs = row_start[i];
    re = row_start[i+1];

    c = num_vertices;
    while(rs < re) {
      if (c != col_index[rs]) {
        // This is an unseen nonzero

        c = col_index[rs];
        ++mRowStart[i];
        ++nnz;
      }

      if (col_instr[rs] >= 0) {
        mAccumulation[col_instr[rs]] = nnz - 1;
      }
      else {
        mAccumulation[-col_instr[rs]] = 1 - nnz;
      }
      
      ++rs;
    }
  }

  nnz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = mRowStart[i];
    mRowStart[i] = nnz;
    nnz += j;
  }
  mRowStart[i] = nnz;

  delete [] col_instr;

  // Fill in the compacted hessian matrix

  mColIndex = new size_t[nnz];

  for (i = 0; i < num_vertices; ++i) {
    rs = row_start[i];
    re = row_start[i+1];

    c = num_vertices;
    while(rs < re) {
      if (c != col_index[rs]) {
        // This is an unseen nonzero

        c = col_index[rs];
        mColIndex[mRowStart[i]] = c;
        mRowStart[i]++;
      }
      ++rs;
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    mRowStart[i+1] = mRowStart[i];
  }
  mRowStart[1] = mRowStart[0];
  mRowStart[0] = 0;
  
  delete [] row_start;
  delete [] col_index;

  mEntries = new Matrix3D[nnz]; // On Solaris, no initializer allowed for new of an array 
  for (i=0;i<nnz;++i) mEntries[i] = 0.; // so we initialize all entries manually. 

  origin_pd = &pd;

  FUNCTION_TIMER_END();
  return;
}


#undef __FUNC__
#define __FUNC__ "MsqHessian::get_diagonal_blocks"
/*! \param diag is an STL vector of size MsqHessian::size() . */
void MsqHessian::get_diagonal_blocks(std::vector<Matrix3D> &diag, MsqError &err)
{
  // make sure we have enough memory, so that no reallocation is needed later.
  if (diag.size() != size()) {
    diag.reserve(size());
  }

  for (size_t i=0; i<size(); ++i) {
    diag[i] = mEntries[mRowStart[i]];
  }
}


#undef __FUNC__
#define __FUNC__ "MsqHessian::compute_preconditionner"
/*! compute a preconditionner used in the preconditionned conjugate gradient
  algebraic solver. In fact, this computes \f$ M^{-1} \f$ .
*/
void MsqHessian::compute_preconditionner(MsqError &err)
{
  // reallocates arrays if size of the Hessian has changed too much.
  if (mSize > precondArraySize || mSize < precondArraySize/10 ) {
    delete[] mPreconditionner;
    mPreconditionner = new Matrix3D[mSize];
  }

  Matrix3D* diag_block;
  double sum, inv_sum;
  size_t m;
  // For each diagonal block, the (inverted) preconditionner is
  // the inverse of the sum of the diagonal entries.
  for (m=0; m<mSize; ++m) {
    diag_block = mEntries + mRowStart[m]; // Gets block at position m,m .

    // find sum, and computes inverse, or 0 if sum = 0 .
    sum = (*diag_block)[0][0] + (*diag_block)[1][1] + (*diag_block)[2][2];
    if (sum != 0.) 
      inv_sum = 1 / sum;
    else
      inv_sum = 0.;
    
    mPreconditionner[m][0][0] = inv_sum;
    mPreconditionner[m][1][1] = inv_sum;
    mPreconditionner[m][2][2] = inv_sum;
  }
}


#undef __FUNC__
#define __FUNC__ "MsqHessian::cg_solver"
/*! uses the preconditionned conjugate gradient algebraic solver
  to find d in \f$ H * d = -g \f$ .
  \param x : the solution, usually the descent direction d.
  \param b : -b will be the right hand side. Usually b is the gradient.
*/
void MsqHessian::cg_solver(Vector3D x[], Vector3D b[], MsqError &err)
{
  FUNCTION_TIMER_START(__FUNC__);
  
  // reallocates arrays if size of the Hessian has changed too much.
  if (mSize > cgArraySizes || mSize < cgArraySizes/10 ) {
    delete[] r;
    delete[] z;
    delete[] p;
    delete[] w;
    r = new Vector3D[mSize];
    z = new Vector3D[mSize];
    p = new Vector3D[mSize];
    w = new Vector3D[mSize];
    cgArraySizes = mSize;
  }

  size_t i;
  double alpha_, alpha, beta; 
  double cg_tol =10e-2; // 10e-2 will give a reasonably good solution (~1%). 
  double norm_g = length(b, mSize);
  double norm_r = norm_g;
  double rzm1; // r^T_{k-1} z_{k-1}
  double rzm2; // r^T_{k-2} z_{k-2}
  this->compute_preconditionner(err); MSQ_CHKERR(err); // get M^{-1} for diagonal blocks

  for (i=0; i<mSize; ++i)  x[i] = 0. ;  
  for (i=0; i<mSize; ++i)  r[i] = -b[i] ;  // r = -b because x_0 = 0 and we solve H*x = -b
  norm_g *= cg_tol;

  this->apply_preconditionner(z, r, err); // solve Mz = r (computes z = M^-1 r)
  for (i=0; i<mSize; ++i)  p[i] = z[i] ; // p_1 = z_0  
  rzm1 = inner(z,r,mSize); // inner product r_{k-1}^T z_{k-1} 
    
  size_t cg_iter = 0;
  while ((norm_r > norm_g) && (cg_iter < maxCGiter)) {
    ++cg_iter;
      
    axpy(w, mSize, *this, p, mSize, 0,0,err); // w = A * p_k
      
    alpha_ = inner(p,w,mSize); // alpha_ = p_k^T A p_k
    if (alpha_ <= 0.0) {
      printf("Direction of Negative Curvature\n");
      break; // Newton goes on with this direction of negative curvature 
    }
      
    alpha = rzm1 / alpha_;
      
    for (i=0; i<mSize; ++i)  x[i] += alpha*p[i]; // x_{k+1} = x_k + alpha_{k+1} p_{k+1} 
    for (i=0; i<mSize; ++i)  r[i] -= alpha*w[i]; // r_{k+1} = r_k - alpha_{k+1} A p_{k+1} 
    norm_r = length(r, mSize);
      
    this->apply_preconditionner(z, r, err); // solve Mz = r (computes z = M^-1 r)
      
    rzm2 = rzm1;
    rzm1 = inner(z,r,mSize); // inner product r_{k-1}^T z_{k-1} 
    beta = rzm1 / rzm2;
    for (i=0; i<mSize; ++i)  p[i] = z[i] + beta*p[i]; // p_k = z_{k-1} + Beta_k * p_{k-1}
  }

  FUNCTION_TIMER_END();
}
