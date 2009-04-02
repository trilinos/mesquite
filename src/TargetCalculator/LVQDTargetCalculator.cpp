/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file LVQDTargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LVQDTargetCalculator.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

LVQDTargetCalculator::LVQDTargetCalculator( 
                        TargetCalculator* lambda_source,
                        TargetCalculator* V_source,
                        TargetCalculator* Q_source,
                        TargetCalculator* D_source )
  : lambdaGuide( lambda_source ),
    vGuide( V_source ),
    qGuide( Q_source ),
    dGuide( D_source )
{ }

LVQDTargetCalculator::~LVQDTargetCalculator() {}

double LVQDTargetCalculator::calc_lambda_3D( const MsqMatrix<3,3>& M )
{
  return Mesquite::cbrt( fabs( det( M ) ) );
}

MsqMatrix<3,3> LVQDTargetCalculator::calc_V_3D( const MsqMatrix<3,3>& M )
{
  MsqMatrix<3,1> c[3] = { M.column(0), M.column(1), M.column(2) };
  const double dotc0c1 = c[0] % c[1];
  MsqMatrix<3,1> c0xc1 = c[0] * c[1];
  const double inv_lenp = 1.0/length(c0xc1);
  const double inv_len0 = 1.0/length(c[0]);
  c[1] *= (c[0] % c[0]);
  c[1] -= dotc0c1*c[0];
  c[1] *= inv_len0 * inv_lenp;
  c[0] *= inv_len0;
  c[2] = inv_lenp * c0xc1;
  return MsqMatrix<3,3>(c);
}

MsqMatrix<3,3> LVQDTargetCalculator::calc_Q_3D( const MsqMatrix<3,3>& M )
{
  const MsqMatrix<3,1> c[3] = { M.column(0),  M.column(1),  M.column(2) };
  const double len[3] = { length(c[0]), length(c[1]), length(c[2]) };
  const MsqMatrix<3,1> c0xc1 = c[0] * c[1], c0xc2 = c[0] * c[2];
  const double len01 = length(c0xc1);
  const double d = det(M);
  const double f = Mesquite::cbrt( len[0] * len[1] * len[2] / fabs(d) );
  
  MsqMatrix<3,3> Q;
  Q(0,0) = f;
  Q(0,1) = f * (c[0] % c[1]) / (len[0] * len[1]);
  Q(0,2) = f * (c[0] % c[2]) / (len[0] * len[2]);
  Q(1,0) = 0.0;
  Q(1,1) = f * len01 / (len[0] * len[1]);
  Q(1,2) = f * (c0xc1 % c0xc2) / (len01 * len[0] * len[2]);
  Q(2,0) = Q(2,1) = 0.0;
  Q(2,2) = f * d / (len01 * len[2]);
  return Q;
}


MsqMatrix<3,3> LVQDTargetCalculator::calc_delta_3D( const MsqMatrix<3,3>& M )
{
  const MsqMatrix<3,1> c[3] = { M.column(0),  M.column(1),  M.column(2) };
  const double len[3] = { length(c[0]), length(c[1]), length(c[2]) };
  const double f = 1.0 / Mesquite::cbrt( len[0] * len[1] * len[2] );
  MsqMatrix<3,3> D(0.0);
  D(0,0) = f * len[0];
  D(1,1) = f * len[1];
  D(2,2) = f * len[2];
  return D;
}

bool LVQDTargetCalculator::get_3D_target( PatchData& pd, 
                                          size_t element,
                                          Sample sample,
                                          MsqMatrix<3,3>& W_out,
                                          MsqError& err )
{
  lambdaGuide->get_3D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const double lambda = calc_lambda_3D( W_out );
  
  vGuide->get_3D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const MsqMatrix<3,3> V = calc_V_3D( W_out );
  
  qGuide->get_3D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const MsqMatrix<3,3> Q = calc_Q_3D( W_out );
  
  dGuide->get_3D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const MsqMatrix<3,3> D = calc_delta_3D( W_out );
  
  W_out = lambda * V * Q * D;
  return true;
}

double LVQDTargetCalculator::calc_lambda_2D( const MsqMatrix<3,2>& M )
{
  const double d = det( transpose(M)*M );
  return sqrt( sqrt( d ) );
}

MsqMatrix<3,2> LVQDTargetCalculator::calc_V_2D( const MsqMatrix<3,2>& M )
{
  MsqMatrix<3,1> c[] = { M.column(0), M.column(1) };
  const double lambda2 = sqrt( det( transpose(M) * M ) );
  const double len0 = length(c[0]);
  const double c0dotc1 = c[0] % c[1];
  c[1] *= c[0]%c[0];
  c[1] -= c0dotc1 * c[0];
  c[1] /= (lambda2 * len0);
  c[0] /= len0;
  return MsqMatrix<3,2>(c);
}

MsqMatrix<2,2> LVQDTargetCalculator::calc_Q_2D( const MsqMatrix<3,2>& M )
{
  const MsqMatrix<3,1> c[] = { M.column(0), M.column(1) };
  const double lambda2 = sqrt( det( transpose(M) * M ) );
  const double lp = sqrt( (c[0] % c[0]) * (c[1] % c[1]) );
  const double q = sqrt( lp / lambda2 );
  MsqMatrix<2,2> Q;
  Q(0,0) = q;
  Q(1,0) = 0.0;
  Q(0,1) = c[0] % c[1] * q/lp;
  Q(1,1) = 1.0/q;
  return Q;
}

MsqMatrix<2,2> LVQDTargetCalculator::calc_delta_2D( const MsqMatrix<3,2>& M )
{
  const MsqMatrix<3,1> c1 = M.column(0);
  const MsqMatrix<3,1> c2 = M.column(1);
  const double r1 = sqrt(length(c1)), r2 = sqrt(length(c2));
  MsqMatrix<2,2> D;
  D(0,0) = r1 / r2;
  D(0,1) = D(1,0) = 0.0;
  D(1,1) = r2 / r1;
  return D;
}


bool LVQDTargetCalculator::get_2D_target( PatchData& pd,
                                          size_t element,
                                          Sample sample,
                                          MsqMatrix<3,2>& W_out,
                                          MsqError& err )
{
  lambdaGuide->get_2D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const double lambda = calc_lambda_2D( W_out );
  
  vGuide->get_2D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const MsqMatrix<3,2> V = calc_V_2D( W_out );
  
  qGuide->get_2D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const MsqMatrix<2,2> Q = calc_Q_2D( W_out );
  
  dGuide->get_2D_target( pd, element, sample, W_out, err ); MSQ_ERRZERO(err);
  const MsqMatrix<2,2> D = calc_delta_2D( W_out );
  
  W_out = lambda * V * Q * D;
  return true;
}


} // namespace Mesquite
