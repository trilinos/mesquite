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


/** \file TargetSkew.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_SKEW_HPP
#define MSQ_TARGET_SKEW_HPP

#include "Mesquite.hpp"
#include "MsqMatrix.hpp"
#include "Sample.hpp"

namespace MESQUITE_NS {

class PatchData;
class MsqError;

/**\brief Class to acquire skew component of target matrices */
class TargetSkew
{
  public:
    virtual
    bool get_skew_2D( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<2,2>& M_out,
                      MsqError& err ) = 0;

    virtual
    bool get_skew_3D( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,3>& M_out,
                      MsqError& err ) = 0;
                        
    
      /**\brief Factor skew component from a 3x3 target matrix */
    inline static 
    bool factor_skew( const MsqMatrix<3,3>& M, MsqMatrix<3,3>& Q_out );
    
      /**\brief Factor skew component from a 3x2 target matrix */
    inline static 
    bool factor_skew( const MsqMatrix<3,2>& M, MsqMatrix<2,2>& Q_out  );
    
      /**\brief Calculate skew matrix from angle */
    inline static
    MsqMatrix<2,2> calc_skew( double phi );
};

inline 
bool TargetSkew::factor_skew( const MsqMatrix<3,3>& M, MsqMatrix<3,3>& Q )
{
  const MsqMatrix<3,1> c[3] = { M.column(0),  M.column(1),  M.column(2) };
  const double len[3] = { length(c[0]), length(c[1]), length(c[2]) };
  const MsqMatrix<3,1> c0xc1 = c[0] * c[1], c0xc2 = c[0] * c[2];
  const double len01 = length(c0xc1);
  const double d = det(M);
  if (fabs(d) < 1e-100)
    return false;
  const double f = Mesquite::cbrt( len[0] * len[1] * len[2] / fabs(d) );
  
  Q(0,0) = f;
  Q(0,1) = f * (c[0] % c[1]) / (len[0] * len[1]);
  Q(0,2) = f * (c[0] % c[2]) / (len[0] * len[2]);
  Q(1,0) = 0.0;
  Q(1,1) = f * len01 / (len[0] * len[1]);
  Q(1,2) = f * (c0xc1 % c0xc2) / (len01 * len[0] * len[2]);
  Q(2,0) = Q(2,1) = 0.0;
  Q(2,2) = f * d / (len01 * len[2]);
  return true;
}

inline 
bool TargetSkew::factor_skew( const MsqMatrix<3,2>& M, MsqMatrix<2,2>& Q )
{
  const MsqMatrix<3,1> c[] = { M.column(0), M.column(1) };
  const double alpha = sqrt( det( transpose(M) * M ) );
  if (alpha < 1e-100)
    return false;
  const double lp = sqrt( (c[0] % c[0]) * (c[1] % c[1]) );
  const double q = sqrt( lp / alpha );
  Q(0,0) = q;
  Q(1,0) = 0.0;
  Q(0,1) = c[0] % c[1] * q/lp;
  Q(1,1) = 1.0/q;
  return true;
}

inline
MsqMatrix<2,2> TargetSkew::calc_skew( double phi )
{
  double sin_phi = sin(phi);
  double f = 1/sqrt(sin_phi);
  MsqMatrix<2,2> result;
  result(0,0) = f;
  result(0,1) = f * cos(phi);
  result(1,0) = 0.0;
  result(1,1) = f * sin_phi;
  return result;
}


} // namespace MESQUITE_NS

#endif
