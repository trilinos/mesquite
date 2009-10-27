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


/** \file TargetAspect.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_ASPECT_HPP
#define MSQ_TARGET_ASPECT_HPP

#include "Mesquite.hpp"
#include "MsqMatrix.hpp"
#include "Sample.hpp"

namespace MESQUITE_NS {

class PatchData;
class MsqError;

/**\brief Class to acquire aspect component of target matrices */
class TargetAspect
{
  public:
    virtual
    bool get_aspect_2D( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<2,2>& M_out,
                        MsqError& err ) = 0;

    virtual
    bool get_aspect_3D( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,3>& M_out,
                        MsqError& err ) = 0;
                        
    
      /**\brief Factor aspect component from a 3x3 target matrix */
    inline static 
    bool factor_aspect( const MsqMatrix<3,3>& M, MsqMatrix<3,3>& D_out );
    
      /**\brief Factor aspect component from a 3x2 target matrix */
    inline static 
    bool factor_aspect( const MsqMatrix<3,2>& M, MsqMatrix<2,2>& D_out );
    
      /**\brief Calculate aspect matrix from vectors */
    inline static
    MsqMatrix<2,2> calc_aspect( MsqVector<2> r );
    
      /**\brief Calculate aspect matrix from vectors */
    inline static
    MsqMatrix<3,3> calc_aspect( MsqVector<3> r );
};

inline 
bool TargetAspect::factor_aspect( const MsqMatrix<3,3>& M, MsqMatrix<3,3>& D )
{
  const MsqMatrix<3,1> c[3] = { M.column(0),  M.column(1),  M.column(2) };
  const double len[3] = { length(c[0]), length(c[1]), length(c[2]) };
  double f;
  if (!divide( 1.0, Mesquite::cbrt( len[0] * len[1] * len[2] ), f ))
    return false;
  D(0,1) = D(0,2) = D(1,0) = D(1,2) = D(2,0) = D(2,1) = 0;
  D(0,0) = f * len[0];
  D(1,1) = f * len[1];
  D(2,2) = f * len[2];
  return true;
}

inline 
bool TargetAspect::factor_aspect( const MsqMatrix<3,2>& M, MsqMatrix<2,2>& D )
{
  const MsqMatrix<3,1> c1 = M.column(0);
  const MsqMatrix<3,1> c2 = M.column(1);
  const double r1 = sqrt(length(c1)), r2 = sqrt(length(c2));
  D(0,1) = D(1,0) = 0.0;
  return divide( r1, r2, D(0,0) ) && divide( r2, r1, D(1,1) );
}

inline
MsqMatrix<2,2> TargetAspect::calc_aspect( MsqVector<2> r )
{
  double root_rho = sqrt(r[0]/r[1]);
  MsqMatrix<2,2> D(0.0);
  D(0,0) = root_rho;
  D(1,1) = 1.0/root_rho;
  return D;
}

inline
MsqMatrix<3,3> TargetAspect::calc_aspect( MsqVector<3> r )
{
  MsqMatrix<3,3> D(0.0);
  D(0,0) = Mesquite::cbrt( r[0]*r[0]/(r[1]*r[2]) );
  D(1,1) = Mesquite::cbrt( r[1]*r[1]/(r[0]*r[2]) );
  D(2,2) = Mesquite::cbrt( r[2]*r[2]/(r[0]*r[1]) );
  return D;
}


} // namespace MESQUITE_NS

#endif
