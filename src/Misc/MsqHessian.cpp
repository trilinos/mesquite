// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 by Thomas Leurent
//  LAST-MOD: 14-Jan-03 at 13:35:26 by Thomas Leurent
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

/*! \brief creates a sparse structure for a Hessian, based on the
  connectivity information contained in the PatchData. */

void MsqHessian::initialize(PatchData &pd, MsqError &err)
{
  MsqMeshEntity* element_array;
  size_t num_vertices = pd.num_vertices();
  if (num_vertices == 0)
    err.set_msg("No vertices in PatchData");
  
  size_t num_elements = pd.num_elements();

  std::vector<size_t> vtx_list;

  size_t e, r, rs, re, c, cs, ce, nz, nnz, nve, i, j;
  element_array = pd.get_element_array(err); MSQ_CHKERR(err);

  // Calculate the offsets for a CSC representation of the accumulation
  // pattern.

  size_t* colStart = new size_t[num_vertices + 1];
    
  for (i = 0; i < num_vertices; ++i) {
    colStart[i] = 0;
  }

  for (e = 0; e < num_elements; ++e) {
    nve = element_array[e].vertex_count();
    element_array[e].get_vertex_indices(vtx_list);

    for (i = 0; i < nve; ++i) {
      colStart[vtx_list[i]] += nve;
    }
  }

  nz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = colStart[i];
    colStart[i] = nz;
    nz += j;
  }
  colStart[i] = nz;

  // Finished putting matrix into CSC representation

  size_t* rowInstr = new size_t[nz];
  size_t* rowIndex = new size_t[nz];

  nz = 0;
  for (e = 0; e < num_elements; ++e) {
    nve = element_array[e].vertex_count();
    element_array[e].get_vertex_indices(vtx_list);

    for (i = 0; i < nve; ++i) {
      r = vtx_list[i];

      for (j = 0 ;j < nve; ++j) {
        c = vtx_list[j];

        rowIndex[colStart[c]] = r;
        rowInstr[colStart[c]] = nz;

        ++colStart[c];
        ++nz;
      }
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    colStart[i+1] = colStart[i];
  }
  colStart[1] = colStart[0];
  colStart[0] = 0;

  // Convert CSC to CSR
  // First calculate the offsets in the row

  size_t* rowStart = new size_t[num_vertices + 1];
    
  for (i = 0; i < num_vertices; ++i) {
    rowStart[i] = 0;
  }

  for (i = 0; i < nz; ++i) {
    ++rowStart[rowIndex[i]];
  }

  nz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = rowStart[i];
    rowStart[i] = nz;
    nz += j;
  }
  rowStart[i] = nz;
    
  // Now calculate the pattern

  size_t* colIndex = new size_t[nz];
  size_t* colInstr = new size_t[nz];

  for (i = 0; i < num_vertices; ++i) {
    cs = colStart[i];
    ce = colStart[i+1];

    while(cs < ce) {
      r = rowIndex[cs];

      colIndex[rowStart[r]] = i;
      colInstr[rowStart[r]] = rowInstr[cs];

      ++rowStart[r];
      ++cs;
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    rowStart[i+1] = rowStart[i];
  }
  rowStart[1] = rowStart[0];
  rowStart[0] = 0;

  delete[] rowIndex;

  // Now that the matrix is CSR
  // Column indices for each row are sorted

  // Compaction -- count the number of nonzeros
  mRowStart = colStart;   // don't need to reallocate
  mColInstr = rowInstr;   // don;t need to reallocate

  for (i = 0; i <= num_vertices; ++i) {
    mRowStart[i] = 0;
  }

  nnz = 0;
  for (i = 0; i < num_vertices; ++i) {
    rs = rowStart[i];
    re = rowStart[i+1];

    c = num_vertices;
    while(rs < re) {
      if (c != colIndex[rs]) {
        // This is an unseen nonzero

        c = colIndex[rs];
        ++mRowStart[i];
        ++nnz;
      }

      mColInstr[colInstr[rs]] = nnz - 1;
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

  delete [] colInstr;

  // Fill in the compacted hessian matrix

  mColIndex = new size_t[nnz];

  for (i = 0; i < num_vertices; ++i) {
    rs = rowStart[i];
    re = rowStart[i+1];

    c = num_vertices;
    while(rs < re) {
      if (c != colIndex[rs]) {
        // This is an unseen nonzero

        c = colIndex[rs];
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
  
  delete [] rowStart;
  delete [] colIndex;

  mEntries = new Matrix3D[nnz];
  return;
}

