/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TargetOrientation.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_ORIENTATION_HPP
#define MSQ_TARGET_ORIENTATION_HPP

#include "Mesquite.hpp"
#include "MsqMatrix.hpp"
#include "Sample.hpp"

namespace MESQUITE_NS {

class PatchData;
class MsqError;

/**\brief Interface for aquiring orientation component of
 *        target matrices
 */
class TargetOrientation
{
  public:
    virtual
    bool get_orient_2D( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,2>& M_out,
                        MsqError& err ) = 0;

    virtual
    bool get_orient_3D( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,3>& M_out,
                        MsqError& err ) = 0;
                        
    
      /**\brief Factor orientation component from a 3x3 target matrix */
    inline static 
    bool factor_orientation( const MsqMatrix<3,3>& M, MsqMatrix<3,3>& V_out );
    
      /**\brief Factor orientation component from a 3x2 target matrix */
    inline static 
    bool factor_orientation( const MsqMatrix<3,2>& M, MsqMatrix<3,2>& V_out );
    
      /**\brief Calculate orientation matrix from vectors */
    inline static
    MsqMatrix<3,3> calc_orientation( const MsqVector<3>& u,
                                     const MsqVector<3>& v );
};

inline 
bool TargetOrientation::factor_orientation( const MsqMatrix<3,3>& M, MsqMatrix<3,3>& V )
{
  MsqMatrix<3,1> c[3] = { M.column(0), M.column(1), M.column(2) };
  const double dotc0c1 = c[0] % c[1];
  MsqMatrix<3,1> c0xc1 = c[0] * c[1];
  double inv_lenp, inv_len0;
  if (!divide(1.0, length(c0xc1), inv_lenp) ||
      !divide(1.0, length(c[0]),  inv_len0) )
    return false;
  c[1] *= (c[0] % c[0]);
  c[1] -= dotc0c1*c[0];
  c[1] *= inv_len0 * inv_lenp;
  c[0] *= inv_len0;
  c[2] = inv_lenp * c0xc1;
  V = MsqMatrix<3,3>(c);
  return true;
}

inline 
bool TargetOrientation::factor_orientation( const MsqMatrix<3,2>& M, MsqMatrix<3,2>& V )
{
  MsqMatrix<3,1> c[] = { M.column(0), M.column(1) };
  const double lambda2 = sqrt( det( transpose(M) * M ) );
  const double len0 = length(c[0]);
  const double c0dotc1 = c[0] % c[1];
  c[1] *= c[0]%c[0];
  c[1] -= c0dotc1 * c[0];
  double f;
  if (!divide(1.0, lambda2*len0, f))
    return false;
  c[1] *= f;
  c[0] /= len0;
  V = MsqMatrix<3,2>(c);
  return true;
}

inline
MsqMatrix<3,3> TargetOrientation::calc_orientation( const MsqVector<3>& u,
                                                    const MsqVector<3>& v )
{
  double u_len_sqr = u % u;
  double u_len = sqrt(u_len_sqr);
  double inv_u_len = 1.0/u_len;
  MsqVector<3> cross = u * v;
  double inv_cross_len = 1.0/sqrt(cross % cross);
  
  MsqMatrix<3,1> cols[3];
  cols[0] = u;
  cols[0] *= inv_u_len;
  cols[1] = v;
  cols[1] *= u_len_sqr;
  cols[1] -= (u % v) * u;
  cols[1] *= inv_u_len * inv_cross_len;
  cols[2] = cross;
  cols[2] *= inv_cross_len;
  return MsqMatrix<3,3>(cols);
}


} // namespace MESQUITE_NS

#endif
