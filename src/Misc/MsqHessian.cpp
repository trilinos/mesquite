// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 by Thomas Leurent
//  LAST-MOD: 22-Jan-03 at 14:28:05 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqHessian.cpp

describe MsqHessian.cpp here

*/


#include "MsqHessian.hpp"

using namespace Mesquite;

using std::cout;
using std::endl;

#undef __FUNC__
#define __FUNC__ "MsqHessian::initialize"
/*! \brief creates a sparse structure for a Hessian, based on the
  connectivity information contained in the PatchData.
  Only the upper triangular part of the Hessian is stored. */
void MsqHessian::initialize(PatchData &pd, MsqError &err)
{
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
    
  for (i = 0; i < num_vertices; ++i) {
    col_start[i] = 0;
  }

  for (e = 0; e < num_elements; ++e) {
    nve = element_array[e].vertex_count();
    element_array[e].get_vertex_indices(vtx_list);

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

  cout << "col_start: ";
  for (int t=0; t<num_vertices+1; ++t)
    cout << col_start[t] << " ";
  cout << endl;
  
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

  cout << "col_start: ";
  for (int t=0; t<num_vertices+1; ++t)
    cout << col_start[t] << " ";
  cout << endl;
  cout << "row_index: ";
  for (int t=0; t<nz; ++t)
    cout << row_index[t] << " ";
  cout << endl;
  cout << "row_instr: ";
  for (int t=0; t<nz; ++t)
    cout << row_instr[t] << " ";
  cout << endl;
  
  
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
  mAccumulation = row_instr;   // don;t need to reallocate

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

  mEntries = new Matrix3D[nnz];
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

  for (int i=0; i<size(); ++i) {
    diag[i] = mEntries[mRowStart[i]];
  }
}


#undef __FUNC__
#define __FUNC__ "axpy"
/*!
  \param res: array of Vector3D in which the result is stored.
  \param x: vector multiplied by the Hessian.
  \param y: vector added to the Hessian vector product.
*/
void axpy(Vector3D res[], int size_r,
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
  int* col = H.mColIndex;
  
  for(i=0; i<nn; ++i)
    res[i] = y[i];
   
  for (i = 0; i < nn; ++i) {
    rl = H.mRowStart[i+1] - H.mRowStart[i]; 
    e = H.mRowStart[i];
    if (rl) {
       
      // Diagonal entry
      res[i] += H.mEntries[e++]*x[i];
      assert(*col == i); // dbg

      //Non-diagonal entries
      for (j = 1; j < rl; ++j) {
        c = *(++col);
         
        res[i] += H.mEntries[e] * x[c];
        res[c] += transpose(H.mEntries[e]) * x[i];
      }
      ++col;
    }
  }
}
