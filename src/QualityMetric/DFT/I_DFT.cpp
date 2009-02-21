/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   I_DFT.cpp
  \brief  

  \author Thomas Leurent

  \date   2004-04-12
*/

#include "I_DFT.hpp"
#include "I_DFTFamilyFunctions.hpp"
#include "TargetMatrix.hpp"

using namespace Mesquite;
   

bool I_DFT::evaluate(PatchData& pd,
			     size_t idx,
			     double& m, 
			     MsqError &err)
{
  // Only works with the weighted average

  MsqMeshEntity* e = &pd.element_by_index( idx );
  const MsqVertex *vertices = pd.get_vertex_array(err); MSQ_ERRZERO(err);

  EntityTopology topo = e->get_element_type();

  //const size_t nv = e->vertex_count();
  const size_t *v_i = e->get_vertex_index_array();

  get_W_matrices( idx, pd, W, 8, mCk, err );
  MSQ_ERRZERO(err);

  // Initialize constants for the metric
  const double delta = useBarrierDelta ? pd.get_barrier_delta(err) :
    (mGamma.value() ? 0 : 1);
  MSQ_ERRZERO(err);
  
  const int triInd[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 0, 3, 2},
                            {2, 3, 0, 1}, {3, 2, 1, 0}};
  const int pyrInd[4][4] = {{0, 1, 3, 4}, {1, 2, 0, 4},
			    {2, 3, 1, 4}, {3, 0, 2, 4}};
  const int priInd[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
			    {2, 0, 1, 5}, {3, 5, 4, 0},
			    {4, 3, 5, 1}, {5, 4, 3, 2}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j;

  m = 0.0;
  switch(topo) {
  case TRIANGLE:
    //assert(3 == nv);

    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[triInd[i][j]]];
      }
      
      mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_2(mMetric, mCoords, mNormals[i], invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      if (!mValid) return false;
      m += mCk[i] * mMetric;
      
    }

    m *= MSQ_ONE_THIRD;
    break;

  case QUADRILATERAL:
    //assert(4 == nv);

    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_2(mMetric, mCoords, mNormals[i], invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      if (!mValid) return false;
      m += mCk[i] * mMetric;
    }

    m *= 0.25;
    break;

  case TETRAHEDRON:
    //assert(4 == nv);

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      
      if (!mValid) return false;
      m += mCk[i] * mMetric;
    }

    m *= 0.25;
    break;

  case PYRAMID:
    //assert(5 == nv);

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[pyrInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      
      if (!mValid) return false;
      m += mCk[i] * mMetric;
    }

    m *= 0.25;
    break;

  case PRISM:
    //assert(6 == nv);

    for (i = 0; i < 6; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[priInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      
      if (!mValid) return false;
      m += mCk[i] * mMetric;
    }

    m *= 1.0 / 6.0;
    break;

  case HEXAHEDRON:
    //assert(8 == nv);

    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      
      if (!mValid) return false;
      m += mCk[i] * mMetric;
    }

    m *= 0.125;
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.",MsqError::UNSUPPORTED_ELEMENT);
    return false;
  }

  return true;
}

bool I_DFT::evaluate_with_gradient( PatchData& pd,
                                    size_t idx,
                                    double& m,
                                    msq_std::vector<size_t>& fv,
                                    msq_std::vector<Vector3D>& g,
                                    MsqError& err )
{
  // Only works with the weighted average

  MsqMeshEntity* e = &pd.element_by_index( idx );
  const MsqVertex *vertices = pd.get_vertex_array(err); MSQ_ERRZERO(err);

  EntityTopology topo = e->get_element_type();

  //const size_t nv = e->vertex_count();
  const size_t *v_i = e->get_vertex_index_array();

  get_W_matrices( idx, pd, W, 8, mCk, err );
  MSQ_ERRZERO(err);
  
  const int nv = e->vertex_count();
  unsigned mFree[8];
  fv.clear();
  int nfv = 0;
  for (int ii = 0; ii < nv; ++ii) 
    if (v_i[ii] < pd.num_free_vertices()) {
      fv.push_back(v_i[ii]);
      mFree[nfv++] = ii;
    }
  g.resize(nfv);

  // Initialize constants for the metric
  const double delta = useBarrierDelta ? pd.get_barrier_delta(err) :
    (mGamma.value() ? 0 : 1);
  MSQ_ERRZERO(err);
  
  const int triInd[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 0, 3, 2},
                            {2, 3, 0, 1}, {3, 2, 1, 0}};
  const int pyrInd[4][4] = {{0, 1, 3, 4}, {1, 2, 0, 4},
			    {2, 3, 1, 4}, {3, 0, 2, 4}};
  const int priInd[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
			    {2, 0, 1, 5}, {3, 5, 4, 0},
			    {4, 3, 5, 1}, {5, 4, 3, 2}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j, mVert;
  
  

  m = 0.0;
  switch(topo) {
  case TRIANGLE:
    //assert(3 == nv);
    
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.
    
    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.
      g[0] = 0.0;
      for (i = 0; i < 3; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[triInd[i][j]]];
	  if (v_i[triInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}

	if (mVert >= 0) {
	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_2_v0(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_2_v1(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_2_v2(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  g[0] += mCk[i] * mGrads[0];
	}
	else {
	  // For triangles, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= MSQ_ONE_THIRD;
      g[0] *= MSQ_ONE_THIRD;
    }
    else {
      for (i = 0; i < 3; ++i) {
	mAccGrads[i] = 0.0;
      }
    
      for (i = 0; i < 3; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[triInd[i][j]]];
	}
      
	mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
      
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_2(mMetric, mGrads, mCoords, mNormals[i], invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);

	if (!mValid) return false;
	m += mCk[i] * mMetric;
	for (j = 0; j < 3; ++j) {
	  mAccGrads[triInd[i][j]] += mCk[i] * mGrads[j];
	}
      }

      m *= MSQ_ONE_THIRD;
      for (i = 0; i < 3; ++i) {
	mAccGrads[i] *= MSQ_ONE_THIRD;
      }

      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free 
      // vertices, in the order of fv[].
      for (i = 0; i < nfv; ++i)
        g[i] = mAccGrads[mFree[i]];
    }
    break;

  case QUADRILATERAL:
    //assert(4 == nv);
    
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}

	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_2_v0(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_2_v1(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_2_v2(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  g[0] += mCk[i] * mGrads[0];
	}
	else {
	  // For quadrilaterals, the free vertex only appears in three 
	  // elements.  Therefore, there these accumulations are needed 
	  // to get the true local objective function.  Note: this code 
          // can be commented out for local codes to improve performance 
          // because you are unable to change the contributions from the 
	  // elements where the free vertex does not appear.  (If the 
	  // weight matrices change, then the code needs to be modified.)
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 0.25;
      g[0] *= 0.25;
    }
    else {
      for (i = 0; i < 4; ++i) {
	mAccGrads[i] = 0.0;
      }

      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}

	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_2(mMetric, mGrads, mCoords, mNormals[i], invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);

	if (!mValid) return false;
	m += mCk[i] * mMetric;
	for (j = 0; j < 3; ++j) {
	  mAccGrads[hexInd[i][j]] += mCk[i] * mGrads[j];
	}
      }
    
      m *= 0.25;
      for (i = 0; i < 4; ++i) {
	mAccGrads[i] *= 0.25;
      }

      // This is not very efficient, but is one way to select correct gradients
      // For gradients, info is returned only for free vertices, in the order 
      // of fv[].
      for (i = 0; i < nfv; ++i)
        g[i] = mAccGrads[mFree[i]];
    }
    break;

  case TETRAHEDRON:
    //assert(4 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	  if (v_i[tetInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_3_v0(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_3_v1(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = g_gdft_3_v2(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_3_v3(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  g[0] += mCk[i] * mGrads[0];
	}
	else {
	  // For tetrahedrons, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 0.25;
      g[0] *= 0.25;
    }
    else {
      for (i = 0; i < 4; ++i) {
	mAccGrads[i] = 0.0;
      }
      
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_3(mMetric, mGrads, mCoords, invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);
	
	if (!mValid) return false;
	m += mCk[i] * mMetric;
	
	for (j = 0; j < 4; ++j) {
	  mAccGrads[tetInd[i][j]] += mCk[i] * mGrads[j];
	}
      }
      
      m *= 0.25;
      for (i = 0; i < 4; ++i) {
	mAccGrads[i] *= 0.25;
      }
      
      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free vertices, 
      // in the order of fv[].
      for (i = 0; i < nfv; ++i)
        g[i] = mAccGrads[mFree[i]];
    }
    break;

  case PYRAMID:
    //assert(5 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[pyrInd[i][j]]];
	  if (v_i[pyrInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_3_v0(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_3_v1(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = g_gdft_3_v2(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_3_v3(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  g[0] += mCk[i] * mGrads[0];
	}
	else {
	  // For pyramids, the free vertex does not appear in every element.
	  // Therefore, there these accumulations will not get used.

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 0.25;
      g[0] *= 0.25;
    }
    else {
      for (i = 0; i < 5; ++i) {
	mAccGrads[i] = 0.0;
      }
      
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[pyrInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_3(mMetric, mGrads, mCoords, invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);
	
	if (!mValid) return false;
	m += mCk[i] * mMetric;
	
	for (j = 0; j < 4; ++j) {
	  mAccGrads[pyrInd[i][j]] += mCk[i] * mGrads[j];
	}
      }
      
      m *= 0.25;
      for (i = 0; i < 5; ++i) {
	mAccGrads[i] *= 0.25;
      }
      
      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free vertices, 
      // in the order of fv[].

      for (i = 0; i < nfv; ++i)
        g[i] = mAccGrads[mFree[i]];
    }
    break;

  case PRISM:
    //assert(6 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 6; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[priInd[i][j]]];
	  if (v_i[priInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}

	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_3_v0(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_3_v1(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = g_gdft_3_v2(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_3_v3(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }
	  
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  g[0] += mCk[i] * mGrads[0];
	}
	else {
	  // For prisms, the free vertex only appears in four elements.
	  // Therefore, there these accumulations are needed to get the
	  // true local objective function.  Note: this code can be commented 
	  // out for local codes to improve performance because you are 
	  // unable to change the contributions from the elements where the 
	  // free vertex does not appear.  (If the weight matrices change, 
	  // then the code needs to be modified.)

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }
      
      m *= 1.0 / 6.0;
      g[0] *= 1.0 / 6.0;
    }
    else {
      for (i = 0; i < 6; ++i) {
	mAccGrads[i] = 0.0;
      }
      
      for (i = 0; i < 6; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[priInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_3(mMetric, mGrads, mCoords, invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);
	
	if (!mValid) return false;
	m += mCk[i] * mMetric;
	
	for (j = 0; j < 4; ++j) {
	  mAccGrads[priInd[i][j]] += mCk[i] * mGrads[j];
	}
      }
      
      m *= 1.0 / 6.0;
      for (i = 0; i < 6; ++i) {
	mAccGrads[i] *= 1.0 / 6.0;
      }
      
      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free 
      // vertices, in the order of fv[].
      
      for (i = 0; i < nfv; ++i)
        g[i] = mAccGrads[mFree[i]];
    }
    break;

  case HEXAHEDRON:
    //assert(8 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 8; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}

	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_3_v0(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_3_v1(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = g_gdft_3_v2(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_3_v3(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }
	  
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  g[0] += mCk[i] * mGrads[0];
	}
	else {
	  // For hexahedrons, the free vertex only appears in four elements.
	  // Therefore, there these accumulations are needed to get the
	  // true local objective function.  Note: this code can be commented 
	  // out for local codes to improve performance because you are 
	  // unable to change the contributions from the elements where the 
	  // free vertex does not appear.  (If the weight matrices change, 
	  // then the code needs to be modified.)

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }
      
      m *= 0.125;
      g[0] *= 0.125;
    }
    else {
      for (i = 0; i < 8; ++i) {
	mAccGrads[i] = 0.0;
      }
      
      for (i = 0; i < 8; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_3(mMetric, mGrads, mCoords, invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);
	
	if (!mValid) return false;
	m += mCk[i] * mMetric;
	
	for (j = 0; j < 4; ++j) {
	  mAccGrads[hexInd[i][j]] += mCk[i] * mGrads[j];
	}
      }
      
      m *= 0.125;
      for (i = 0; i < 8; ++i) {
	mAccGrads[i] *= 0.125;
      }
      
      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free 
      // vertices, in the order of fv[].
      
      for (i = 0; i < nfv; ++i)
        g[i] = mAccGrads[mFree[i]];
    }
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.",MsqError::UNSUPPORTED_ELEMENT);
    return false;
  }

  return true;
}

bool I_DFT::evaluate_with_Hessian( PatchData& pd,
                                   size_t idx,
                                   double& m,
                                   msq_std::vector<size_t>& fv,
                                   msq_std::vector<Vector3D>& g,
                                   msq_std::vector<Matrix3D>& h,
                                   MsqError& err )
{
  // Only works with the weighted average

  MsqMeshEntity* e = &pd.element_by_index( idx );
  const MsqVertex *vertices = pd.get_vertex_array(err); MSQ_ERRZERO(err);

  EntityTopology topo = e->get_element_type();

  //const size_t nv = e->vertex_count();
  const size_t *v_i = e->get_vertex_index_array();

  get_W_matrices( idx, pd, W, 8, mCk, err );
  MSQ_ERRZERO(err);
  
  fv.clear();
  uint32_t fixed_bits = fixed_vertex_bitmap( pd, e, fv );
  const int nfv = fv.size();

  const int nv = e->vertex_count();
  g.clear();
  g.resize(nv, Vector3D(0.0));
  h.clear();
  h.resize(nv*(nv+1)/2, Matrix3D(0.0));

  // Initialize constants for the metric
  const double delta = useBarrierDelta ? pd.get_barrier_delta(err) :
    (mGamma.value() ? 0 : 1);  
  MSQ_ERRZERO(err);

  const int triInd[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 0, 3, 2},
                            {2, 3, 0, 1}, {3, 2, 1, 0}};
  const int pyrInd[4][4] = {{0, 1, 3, 4}, {1, 2, 0, 4},
			    {2, 3, 1, 4}, {3, 0, 2, 4}};
  const int priInd[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
			    {2, 0, 1, 5}, {3, 5, 4, 0},
			    {4, 3, 5, 1}, {5, 4, 3, 2}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j, k, l, mVert;
  int      row, col, loc;

  m = 0.0;
  switch(topo) {
  case TRIANGLE:
    //assert(3 == nv);
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 3; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[triInd[i][j]]];
	  if (v_i[triInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_2_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_2_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_2_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  mG += mCk[i] * mGrads[0];
	  mH += mCk[i] * mHessians[0];
	}
	else {
	  // For triangles, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= MSQ_ONE_THIRD;
      mG *= MSQ_ONE_THIRD;
      mH *= MSQ_ONE_THIRD;
      g.resize(1);
      g[0] = mG;
      h.resize(1);
      h[0] = mH;
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 3; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[triInd[i][j]]];
	}

	mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
      
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_2(mMetric, mGrads, mHessians, mCoords, mNormals[i], 
			  invR, mQ, mAlpha, mGamma, delta, mBeta);

	if (!mValid) return false;
	m += mCk[i] * mMetric;

	for (j = 0; j < 3; ++j) {
	  g[triInd[i][j]] += mCk[i] * mGrads[j];
	}

	l = 0;
	for (j = 0; j < 3; ++j) {
	  for (k = j; k < 3; ++k) {
	    row = triInd[i][j];
	    col = triInd[i][k];

	    if (row <= col) {
	      loc = 3*row - (row*(row+1)/2) + col;
	      h[loc] += mCk[i] * mHessians[l];
	    }
	    else {
	      loc = 3*col - (col*(col+1)/2) + row;
	      h[loc] += mCk[i] * transpose(mHessians[l]);
	    }
	    ++l;
	  }
	}
      }

      m *= MSQ_ONE_THIRD;
      for (i = 0; i < 3; ++i) {
	g[i] *= MSQ_ONE_THIRD;
      }

      for (i = 0; i < 6; ++i) {
	h[i] *= MSQ_ONE_THIRD;
      }
      
      if (fixed_bits) {
        remove_fixed_gradients( TRIANGLE, fixed_bits, g );
        remove_fixed_hessians( TRIANGLE, fixed_bits, h );
      }
    }
    break;

  case QUADRILATERAL:
    //assert(4 == nv);
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_2_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_2_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_2_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  mG += mCk[i] * mGrads[0];
	  mH += mCk[i] * mHessians[0];
	}
	else {
	  // For quadrilaterals, the free vertex only appears in three 
	  // elements.  Therefore, there these accumulations are needed 
	  // to get the true local objective function.  Note: this code 
          // can be commented out for local codes to improve performance 
          // because you are unable to change the contributions from the 
	  // elements where the free vertex does not appear.  (If the 
	  // weight matrices change, then the code needs to be modified.)
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 0.25;
      mG *= 0.25;
      mH *= 0.25;
      g.resize(1);
      g[0] = mG;
      h.resize(1);
      h[0] = mH;
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}

	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_2(mMetric, mGrads, mHessians, mCoords, mNormals[i], 
			  invR, mQ, mAlpha, mGamma, delta, mBeta);

	if (!mValid) return false;
	m += mCk[i] * mMetric;

	for (j = 0; j < 3; ++j) {
	  g[hexInd[i][j]] += mCk[i] * mGrads[j];
	}

	l = 0;
	for (j = 0; j < 3; ++j) {
	  for (k = j; k < 3; ++k) {
	    row = hexInd[i][j];
	    col = hexInd[i][k];

	    if (row <= col) {
	      loc = 4*row - (row*(row+1)/2) + col;
	      h[loc] += mCk[i] * mHessians[l];
	    }
	    else {
	      loc = 4*col - (col*(col+1)/2) + row;
	      h[loc] += mCk[i] * transpose(mHessians[l]);
	    }
	    ++l;
	  }
	}
      }

      m *= 0.25;
      for (i = 0; i < 4; ++i) {
	g[i] *= 0.25;
      }

      for (i = 0; i < 10; ++i) {
	h[i] *= 0.25;
      }

      if (fixed_bits) {
        remove_fixed_gradients( QUADRILATERAL, fixed_bits, g );
        remove_fixed_hessians( QUADRILATERAL, fixed_bits, h );
      }
    }
    break;

  case TETRAHEDRON:
    //assert(4 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	  if (v_i[tetInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_3_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_3_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = h_gdft_3_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_3_v3(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  mG += mCk[i] * mGrads[0];
	  mH += mCk[i] * mHessians[0];
	}
	else {
	  // For tetrahedrons, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 0.25;
      mG *= 0.25;
      mH *= 0.25;
      g.resize(1);
      g[0] = mG;
      h.resize(1);
      h[0] = mH;
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_3(mMetric, mGrads, mHessians, mCoords, invR, mQ,
			  mAlpha, mGamma, delta, mBeta);
	
	if (!mValid) return false;
	
	m += mCk[i] * mMetric;
	
	for (j = 0; j < 4; ++j) {
	  g[tetInd[i][j]] += mCk[i] * mGrads[j];
	}
	
	l = 0;
	for (j = 0; j < 4; ++j) {
	  for (k = j; k < 4; ++k) {
	    row = tetInd[i][j];
	    col = tetInd[i][k];
	    
	    if (row <= col) {
	      loc = 4*row - (row*(row+1)/2) + col;
	      h[loc] += mCk[i] * mHessians[l];
	    }
	    else {
	      loc = 4*col - (col*(col+1)/2) + row;
	      h[loc] += mCk[i] * transpose(mHessians[l]);
	    }
	    ++l;
	  }
	}
      }
      
      m *= 0.25;
      for (i = 0; i < 4; ++i) {
	g[i] *= 0.25;
      }
      
      for (i = 0; i < 10; ++i) {
	h[i] *= 0.25;
      }
      
      if (fixed_bits) {
        remove_fixed_gradients( TETRAHEDRON, fixed_bits, g );
        remove_fixed_hessians( TETRAHEDRON, fixed_bits, h );
      }
    }
    break;

  case PYRAMID:
    //assert(5 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[pyrInd[i][j]]];
	  if (v_i[pyrInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_3_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_3_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = h_gdft_3_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_3_v3(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  mG += mCk[i] * mGrads[0];
	  mH += mCk[i] * mHessians[0];
	}
	else {
	  // For pyramids, the free vertex does not appear in every element.
	  // Therefore, there these accumulations will not get used.

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 0.25;
      mG *= 0.25;
      mH *= 0.25;
      g.resize(1);
      g[0] = mG;
      h.resize(1);
      h[0] = mH;
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[pyrInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_3(mMetric, mGrads, mHessians, mCoords, invR, mQ,
			  mAlpha, mGamma, delta, mBeta);
	
	if (!mValid) return false;
	
	m += mCk[i] * mMetric;
	
	for (j = 0; j < 4; ++j) {
	  g[pyrInd[i][j]] += mCk[i] * mGrads[j];
	}
	
	l = 0;
	for (j = 0; j < 4; ++j) {
	  for (k = j; k < 4; ++k) {
	    row = pyrInd[i][j];
	    col = pyrInd[i][k];
	    
	    if (row <= col) {
	      loc = 5*row - (row*(row+1)/2) + col;
	      h[loc] += mCk[i] * mHessians[l];
	    }
	    else {
	      loc = 5*col - (col*(col+1)/2) + row;
	      h[loc] += mCk[i] * transpose(mHessians[l]);
	    }
	    ++l;
	  }
	}
      }
      
      m *= 0.25;
      for (i = 0; i < 5; ++i) {
	g[i] *= 0.25;
      }
      
      for (i = 0; i < 15; ++i) {
	h[i] *= 0.25;
      }
      
      if (fixed_bits) {
        remove_fixed_gradients( PYRAMID, fixed_bits, g );
        remove_fixed_hessians( PYRAMID, fixed_bits, h );
      }
    }
    break;

  case PRISM:
    //assert(6 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 6; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[priInd[i][j]]];
	  if (v_i[priInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_3_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_3_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = h_gdft_3_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_3_v3(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  mG += mCk[i] * mGrads[0];
	  mH += mCk[i] * mHessians[0];
	}
	else {
	  // For prisms, the free vertex only appears in four elements.
	  // Therefore, there these accumulations are needed to get the
	  // true local objective function.  Note: this code can be commented 
	  // out for local codes to improve performance because you are 
	  // unable to change the contributions from the elements where the 
	  // free vertex does not appear.  (If the weight matrices change, 
	  // then the code needs to be modified.)

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 1.0 / 6.0;
      mG *= 1.0 / 6.0;
      mH *= 1.0 / 6.0;
      g.resize(1);
      g[0] = mG;
      h.resize(1);
      h[0] = mH;
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 6; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[priInd[i][j]]];
	}

	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_3(mMetric, mGrads, mHessians, mCoords, invR, mQ,
			  mAlpha, mGamma, delta, mBeta);
      
	if (!mValid) return false;

	m += mCk[i] * mMetric;

	for (j = 0; j < 4; ++j) {
	  g[priInd[i][j]] += mCk[i] * mGrads[j];
	}

	l = 0;
	for (j = 0; j < 4; ++j) {
	  for (k = j; k < 4; ++k) {
	    row = priInd[i][j];
	    col = priInd[i][k];

	    if (row <= col) {
	      loc = 6*row - (row*(row+1)/2) + col;
	      h[loc] += mCk[i] * mHessians[l];
	    }
	    else {
	      loc = 6*col - (col*(col+1)/2) + row;
	      h[loc] += mCk[i] * transpose(mHessians[l]);
	    }
	    ++l;
	  }
	}
      }

      m *= 1.0 / 6.0;
      for (i = 0; i < 6; ++i) {
	g[i] *= 1.0 / 6.0;
      }

      for (i = 0; i < 21; ++i) {
	h[i] *= 1.0 / 6.0;
      }

      if (fixed_bits) {
        remove_fixed_gradients( PRISM, fixed_bits, g );
        remove_fixed_hessians( PRISM, fixed_bits, h );
      }
    }
    break;

  case HEXAHEDRON:
    //assert(8 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 8; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_3_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_3_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = h_gdft_3_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_3_v3(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	  mG += mCk[i] * mGrads[0];
	  mH += mCk[i] * mHessians[0];
	}
	else {
	  // For hexahedrons, the free vertex only appears in four elements.
	  // Therefore, there these accumulations are needed to get the
	  // true local objective function.  Note: this code can be commented 
	  // out for local codes to improve performance because you are 
	  // unable to change the contributions from the elements where the 
	  // free vertex does not appear.  (If the weight matrices change, 
	  // then the code needs to be modified.)

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += mCk[i] * mMetric;
	}
      }

      m *= 0.125;
      mG *= 0.125;
      mH *= 0.125;
      g.resize(1);
      g[0] = mG;
      h.resize(1);
      h[0] = mH;
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 8; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}

	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_3(mMetric, mGrads, mHessians, mCoords, invR, mQ,
			  mAlpha, mGamma, delta, mBeta);
      
	if (!mValid) return false;

	m += mCk[i] * mMetric;

	for (j = 0; j < 4; ++j) {
	  g[hexInd[i][j]] += mCk[i] * mGrads[j];
	}

	l = 0;
	for (j = 0; j < 4; ++j) {
	  for (k = j; k < 4; ++k) {
	    row = hexInd[i][j];
	    col = hexInd[i][k];

	    if (row <= col) {
	      loc = 8*row - (row*(row+1)/2) + col;
	      h[loc] += mCk[i] * mHessians[l];
	    }
	    else {
	      loc = 8*col - (col*(col+1)/2) + row;
	      h[loc] += mCk[i] * transpose(mHessians[l]);
	    }
	    ++l;
	  }
	}
      }

      m *= 0.125;
      for (i = 0; i < 8; ++i) {
	g[i] *= 0.125;
      }

      for (i = 0; i < 36; ++i) {
	h[i] *= 0.125;
      }

      if (fixed_bits) {
        remove_fixed_gradients( HEXAHEDRON, fixed_bits, g );
        remove_fixed_hessians( HEXAHEDRON, fixed_bits, h );
      }
    }
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.",MsqError::UNSUPPORTED_ELEMENT);
    return false;
  }
  return true;
}
