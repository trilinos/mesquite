// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 by Thomas Leurent
//  LAST-MOD: 22-Jan-03 at 17:56:42 by Thomas Leurent
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


MsqHessian::MsqHessian() :
  origin_pd(0), mEntries(0), mRowStart(0), mColIndex(0), 
  mAccumulation(0), mAccumElemStart(0), mSize(0)
{ }


#undef __FUNC__
#define __FUNC__ "MsqHessian::initialize"
/*! \brief creates a sparse structure for a Hessian, based on the
  connectivity information contained in the PatchData.
  Only the upper triangular part of the Hessian is stored. */
void MsqHessian::initialize(PatchData &pd, MsqError &err)
{
  delete[] mEntries;
  delete[] mRowStart;
  delete[] mColIndex;
  delete[] mAccumulation;
  delete[] mAccumElemStart;
  
  size_t num_vertices = pd.num_vertices();
  size_t num_elements = pd.num_elements();
  std::vector<size_t> vtx_list;
  size_t e, r, rs, re, c, cs, ce, nz, nnz, nv, nve, i, j;
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
  nv = 0;
  
  for (i = 0; i < num_vertices; ++i) {
    col_start[i] = 0;
  }

  for (e = 0; e < num_elements; ++e) {
    nve = element_array[e].vertex_count();
    element_array[e].get_vertex_indices(vtx_list);
    nv+=nve;
    mAccumElemStart[e+1] = nv;
    
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
#define __FUNC__ "MsqHessian::accumulate_entries"
/*! \param pd: PatchData in that contains the element which Hessian
           we are accumulating in the Hessian matrix. This must be the same
           PatchData that was used in MsqHessian::initialize().
    \param elem_index: index of the element in the PatchData.
    \param mat3d_array: This is the upper triangular part of the element Hessian for all nodes, including fixed nodes, for which the entries must be null Matrix3Ds.
    \param nb_mat3d. The size of the mat3d_array: (n+1)n/2, where n is
           the number of nodes in the element.
  */
 void MsqHessian::accumulate_entries(PatchData &pd, size_t elem_index,
                         Matrix3D mat3d_array[], size_t nb_mat3d, MsqError &err)
{
  int i;
  
  if (&pd != origin_pd) {
    err.set_msg("Cannot accumulate elements from a different patch. "
                "Use MsqHessian::initialize first.");
    return;
  }

  int nve = pd.get_element_array(err)[elem_index].vertex_count(); MSQ_CHKERR(err);

  assert( nb_mat3d == (nve+1)*nve/2 );

  int e = mAccumElemStart[elem_index];
  for (i=0; i<nb_mat3d; ++i) {
    if (mEntries[e] >= 0)
      mEntries[e] += mat3d_array[i];
    else
      mEntries[e].plus_transpose(mat3d_array[i]);
    ++e;
  }

  assert( e == mAccumElemStart[elem_index+1] );
}
