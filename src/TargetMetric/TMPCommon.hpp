/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TMPCommon.hpp
 *  \brief Common utility stuff for implementing target metrics
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TMP_COMMON_HPP
#define MSQ_TMP_COMMON_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {


#define TMP_T_TEMPL_IMPL_DIM(N,D) \
bool N::evaluate( const MsqMatrix<D,D>& T, double& r, MsqError& ) \
  { return eval( T, r ); } \
bool N::evaluate_with_grad( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqError& ) \
  { return grad( T, r, d1 ); } \
bool N::evaluate_with_Hess( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqMatrix<D,D>* d2, MsqError& ) \
  { return hess( T, r, d1, d2 ); }

#define TMP_T_TEMPL_IMPL_COMMON(N) \
  TMP_T_TEMPL_IMPL_DIM(N,2) \
  TMP_T_TEMPL_IMPL_DIM(N,3) 

#define TMP_T_COMP1_TEMPL_IMPL_DIM(N,D) \
bool N::evaluate( const MsqMatrix<D,D>& T, double& r, MsqError& ) \
  { return eval( mMetric, T, r ); } \
bool N::evaluate_with_grad( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqError& ) \
  { return grad( mMetric, T, r, d1 ); } \
bool N::evaluate_with_Hess( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqMatrix<D,D>* d2, MsqError& ) \
  { return hess( mMetric, T, r, d1, d2 ); }

#define TMP_T_COMP1_TEMPL_IMPL_COMMON(N) \
  TMP_T_COMP1_TEMPL_IMPL_DIM(N,2) \
  TMP_T_COMP1_TEMPL_IMPL_DIM(N,3) 

#define TMP_AW_TEMPL_IMPL_DIM(N,D) \
bool N::evaluate( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqError& ) \
  { return eval( A, W, r ); } \
bool N::evaluate_with_grad( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqMatrix<D,D>& d1, MsqError& ) \
  { return grad( A, W, r, d1 ); } \
bool N::evaluate_with_Hess( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqMatrix<D,D>& d1, MsqMatrix<D,D>* d2, MsqError& ) \
  { return hess( A, W, r, d1, d2 ); }

#define TMP_AW_TEMPL_IMPL_COMMON(N) \
  TMP_AW_TEMPL_IMPL_DIM(N,2) \
  TMP_AW_TEMPL_IMPL_DIM(N,3) 

template <unsigned D> struct DimConst {};
template <> struct DimConst<2> {
  static inline double sqrt() { return MSQ_SQRT_TWO; }
  static inline double inv()  { return 0.5; }
};
template <> struct DimConst<3> {
  static inline double sqrt() { return MSQ_SQRT_THREE; }
  static inline double inv()  { return MSQ_ONE_THIRD; }
};



} // namespace MESQUITE_NS

#endif
