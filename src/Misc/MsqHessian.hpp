// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19
//  LAST-MOD: 10-Jan-03 at 15:30:08 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqHessian.hpp

describe MsqHessian.hpp here

 */


#ifndef MsqHessian_hpp
#define MsqHessian_hpp

#include "Mesquite.hpp"
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
    size_t* mColumnIndex;  //< column indexes of the entries in the row. 
    size_t* mRowStart;	   // size: number of vertices

    
  public:
    void initialize(PatchData &pd, MsqError &err);
    
    
  };

  class MsqHessianIterator
  {
  private:
    Matrix3D* rowStart;
    size_t mPos;
    
  public:
    Matrix3D operator++();
    
  }


  /*! \brief creates a sparse structure for a Hessian, based on the
      connectivity information contained in the PatchData. */
  void MsqHessian::initialize(PatchData &pd, MsqError &err)
  {
    MsqMeshEntity* element_array;
    size_t num_vertices = pd.num_vertices();
    size_t num_elements = pd.num_elements();

    std::vector<size_t> vtx_list;

    size_t e, r, rs, re, c, cs, ce, i, j, nz, nnz;

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

    for (i = num_vertices-1; i >= 0; --i) {
      colStart[i+1] = colStart[i];
    }
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

    nz = 0;
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

    for (i = num_vertices-1; i >= 0; --i) {
      rowStart[i+1] = rowStart[i];
    }
    rowStart[0] = 0;

    delete[] rowIndex;

    // Now have the matrix is CSR
    // Column indieces for each row are sorted

    // Compaction -- count the number of nonzeros
    mRowStart = colStart;   // don't need to reallocate
    mColInstr = rowInstr;   // don;t need to reallocate

    for (i = 0; i <= num_vertices; ++i) {
      mRowStart[i] = 0;
    }

    nz = 0;
    nnz = 0;
    for (i = 0; i < num_vertices; ++i) {
      rs = rowStart[i];
      re = rowStart[i+1];

      c = -1;
      while(rs < re) {
        if (c != colIndex[rs]) {
	  // This is an unseen nonzero

	  c = colIndex[rs];
	  ++rowStart[i];
          ++nnz;
        }

        mColInstr[colInstr[rs]] = nz++;
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

      c = -1;
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

    for (i = num_vertices-1; i >= 0; --i) {
      mRowStart[i+1] = mRowStart[i];
    }

    mRowStart[0] = 0;
    delete [] rowStart;
    delete [] colIndex;

    mEntries = new Matrix3D[nnz];
    return;
  }
  

} // namespace

#endif // MsqHessian_hpp
