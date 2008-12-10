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


/** \file TriLagrangeShape.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TriLagrangeShape.hpp"
#include "MsqError.hpp"
#include <assert.h>

namespace Mesquite {

EntityTopology TriLagrangeShape::element_topology() const
  { return TRIANGLE; }

void TriLagrangeShape::coefficients_at_corner( unsigned corner, 
                                               unsigned nodebits,
                                               msq_std::vector<double>& coeff_out,
                                               MsqError& err) const
{
  if (nodebits >= (1u << 3)) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  coeff_out.clear();
  coeff_out.resize( 6, 0.0 );
  coeff_out[corner] = 1.0;
}

void TriLagrangeShape::coefficients_at_mid_edge( unsigned edge,
                                                 unsigned nodebits,
                                                 msq_std::vector<double>& coeff_out,
                                                 MsqError& err) const
{
  if (nodebits >= (1u << 3)) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  coeff_out.clear();
  coeff_out.resize( 6, 0.0 );
  if (nodebits & (1 << edge)) { // if mid-edge node is present
    coeff_out[3+edge] = 1.0;
  }
  else { // no mid node on edge
    coeff_out[edge] = 0.5;
    coeff_out[(edge+1)%3] = 0.5;
  }
}


void TriLagrangeShape::coefficients_at_mid_face( unsigned ,
                                                unsigned ,
                                                msq_std::vector<double>& ,
                                                MsqError& err) const
{
  MSQ_SETERR(err)("Request for mid-face mapping function value"
                  "for a triangle element", MsqError::UNSUPPORTED_ELEMENT);
}

void TriLagrangeShape::coefficients_at_mid_elem( unsigned nodebits,
                                                 msq_std::vector<double>& coeff_out,
                                                 MsqError& err) const
{
  if (nodebits >= (1u << 3)) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  const double one_third = 1.0/3.0;
  const double two_ninth = 2.0/9.0;
  const double four_ninth = 4.0/9.0;
  
  coeff_out.clear();
  coeff_out.resize( 6, 0.0 );
  
  coeff_out[0] = one_third;
  coeff_out[1] = one_third;
  coeff_out[2] = one_third;
  for (int i = 0; i < 3; ++i) { // for each mid-edge node
    if (nodebits & (1 << i)) { // if node is present
        // coeff for mid-edge node
      coeff_out[i+3]      = four_ninth;
        // adjust coeff for adj corner nodes
      coeff_out[i]       -= two_ninth;
      coeff_out[(i+1)%3] -= two_ninth;
    }
  }
}

static inline void get_linear_derivatives( msq_std::vector<size_t>& vertices,
                                           msq_std::vector<double>& derivs )
{
  vertices.resize(3);
  derivs.resize(6);
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 2;
  derivs[0] =  1.0;
  derivs[1] =  0.0;
  derivs[2] =  0.0;
  derivs[3] =  1.0;
  derivs[4] = -1.0;
  derivs[5] = -1.0;
}  

void TriLagrangeShape::derivatives_at_corner( unsigned corner,
                                              unsigned nodebits,
                                              msq_std::vector<size_t>& vertices,
                                              msq_std::vector<double>& derivs,
                                              MsqError& err) const
{
  if (nodebits >= (1u << 3)) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  get_linear_derivatives( vertices, derivs );
  switch (corner) {
    case 0:
      if (nodebits & 1) {
        vertices.push_back(3);
        derivs.push_back(0.0);
        derivs.push_back(4.0);
        derivs[1] -= 2.0;
        derivs[3] -= 2.0;
      }
      if (nodebits & 4) {
        vertices.push_back(5);
        derivs.push_back(-4.0);
        derivs.push_back(-4.0);
        derivs[0] += 2.0;
        derivs[1] += 2.0;
        derivs[4] += 2.0;
        derivs[5] += 2.0;
      }
      break;
    
    case 1:
      if (nodebits & 1) {
        vertices.push_back(3);
        derivs.push_back(4.0);
        derivs.push_back(0.0);
        derivs[0] -= 2.0;
        derivs[2] -= 2.0;
      }
      if (nodebits & 2) {
        vertices.push_back(4);
        derivs.push_back(-4.0);
        derivs.push_back(-4.0);
        derivs[2] += 2.0;
        derivs[3] += 2.0;
        derivs[4] += 2.0;
        derivs[5] += 2.0;
      }
      break;
    
    case 2:
      if (nodebits & 2) {
        vertices.push_back(4);
        derivs.push_back(0.0);
        derivs.push_back(4.0);
        derivs[3] -= 2.0;
        derivs[5] -= 2.0;
      }
      if (nodebits & 4) {
        vertices.push_back(5);
        derivs.push_back(4.0);
        derivs.push_back(0.0);
        derivs[0] -= 2.0;
        derivs[4] -= 2.0;
      }
      break;
  }
}

static const double edr[3][3] = { { 2.0, 2.0, 0.0 },
                                  {-2.0,-2.0, 0.0 },
                                  {-2.0, 2.0, 0.0 } };
static const double eds[3][3] = { { 2.0, 0.0, 2.0 },
                                  {-2.0, 0.0, 2.0 },
                                  {-2.0, 0.0,-2.0 } };

void TriLagrangeShape::derivatives_at_mid_edge( unsigned edge, 
                                unsigned nodebits,
                                msq_std::vector<size_t>& vertices,
                                msq_std::vector<double>& derivs,
                                MsqError& err ) const
{
  if (nodebits >= (1u << 3)) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
    // The mid-edge behavior is rather strange.
    // A corner vertex contributes to the jacobian
    // at the mid-point of the opposite edge unless
    // one, but *not* both of the adjacent mid-edge
    // nodes is present.
    
    // The value for each corner is incremented when 
    // a mid-side node is present.  If the value is
    // 0 when finished, the corner doesn't contribute.
    // Initialize values to 0 for corners adjacent to
    // edge so they are never zero.
  int corner_count[3] = { 1, 1, 1 };
  corner_count[(edge+2)%3] = -1;
  
    // begin with derivatives for linear tri
  double corner_derivs[3][2] = { { 1.0, 0.0 },
                                 { 0.0, 1.0 },
                                 {-1.0,-1.0 } };
  
  vertices.clear();
  derivs.clear();
  
    // do mid-side nodes first
  for (unsigned i = 0; i < 3; ++i) {
    if (nodebits & (1<<i)) {
      vertices.push_back( i+3 );
      derivs.push_back( edr[i][edge] );
      derivs.push_back( eds[i][edge] );

      int a = (i+1)%3;
      corner_derivs[i][0] -= 0.5 * edr[i][edge];
      corner_derivs[i][1] -= 0.5 * eds[i][edge];
      corner_derivs[a][0] -= 0.5 * edr[i][edge];
      corner_derivs[a][1] -= 0.5 * eds[i][edge];
      ++corner_count[i];
      ++corner_count[a];
    }
  }
  
    // now add corner nodes to list
  for (unsigned i = 0; i < 3; ++i) {
    if (corner_count[i]) {
      vertices.push_back(i);
      derivs.push_back(corner_derivs[i][0]);
      derivs.push_back(corner_derivs[i][1]);
    }
  }
}

void TriLagrangeShape::derivatives_at_mid_face( unsigned , 
                                unsigned ,
                                msq_std::vector<size_t>& ,
                                msq_std::vector<double>& ,
                                MsqError& err ) const
{
  MSQ_SETERR(err)("Request for mid-face mapping function derivative"
                  "for a triangle element", MsqError::UNSUPPORTED_ELEMENT);
}

static const double fdr[] = { 4.0/3.0, -4.0/3.0, 0.0 };
static const double fds[] = { 4.0/3.0, 0.0, -4.0/3.0 };
       
void TriLagrangeShape::derivatives_at_mid_elem( unsigned nodebits,
                                            msq_std::vector<size_t>& vertices,
                                            msq_std::vector<double>& derivs,
                                            MsqError& err ) const
{
  if (nodebits >= (1u << 3)) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  get_linear_derivatives( vertices, derivs );
  for (unsigned i = 0; i < 3; ++i) 
    if (nodebits & (1<<i)) {
      vertices.push_back(i+3);
      derivs.push_back(fdr[i]);
      derivs.push_back(fds[i]);
       
      int a = (i+1)%3;
      derivs[2*i  ] -= 0.5 * fdr[i];
      derivs[2*i+1] -= 0.5 * fds[i];
      derivs[2*a  ] -= 0.5 * fdr[i];
      derivs[2*a+1] -= 0.5 * fds[i];
    }
}
  
} // namespace Mesquite
