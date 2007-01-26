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


/** \file TetLagrangeShape.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TetLagrangeShape.hpp"
#include "MsqError.hpp"
#include <assert.h>

namespace Mesquite {

EntityTopology TetLagrangeShape::element_topology() const
  { return TETRAHEDRON; }

static inline int have_node( unsigned nodebits, unsigned node )
  { return nodebits & (1 << (node-4)); }

void TetLagrangeShape::coefficients_at_corner( unsigned corner, 
                               unsigned nodebits,
                               msq_std::vector<double>& coeff_out,
                               MsqError& err ) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  coeff_out.clear();
  coeff_out.resize( 10, 0.0 );
  coeff_out[corner] = 1.0;
}

void TetLagrangeShape::coefficients_at_mid_edge( unsigned edge, 
                               unsigned nodebits,
                               msq_std::vector<double>& coeff_out,
                               MsqError& err ) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  coeff_out.clear();
  coeff_out.resize( 10, 0.0 );
  if (nodebits & (1 << edge)) { // if mid-edge node is present
    coeff_out[4+edge] = 1.0;
  }
  else { // no mid node on edge
    if (edge < 3) {
      coeff_out[edge] = 0.5;
      coeff_out[(edge+1) % 3] = 0.5;
    }
    else {
      coeff_out[edge-3]= 0.5;
      coeff_out[3] = 0.5;
    }
  }
}

void TetLagrangeShape::coefficients_at_mid_face( unsigned face, 
                               unsigned nodebits,
                               msq_std::vector<double>& coeff_out,
                               MsqError& err ) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  coeff_out.clear();
  coeff_out.resize( 10, 0.0 );
  
  const double one_ninth = 1.0/9.0;
  const double two_ninth = 2.0/9.0;
  const double four_ninth = 4.0/9.0;
  
  if (face < 3) {
    const int next = (face+1)%3;
    coeff_out[face] = -one_ninth;
    coeff_out[next] = -one_ninth;
    coeff_out[3]    = -one_ninth;
    if (nodebits & (1<<face)) {
      coeff_out[4+face] = four_ninth;
    }
    else {
      coeff_out[face] += two_ninth;
      coeff_out[next] += two_ninth;
    }
    if (nodebits & (1<<(3+next))) {
      coeff_out[7+next] = four_ninth;
    }
    else {
      coeff_out[face] += two_ninth;
      coeff_out[3]    += two_ninth;
    }
    if (nodebits & (1<<(3+face))) {
      coeff_out[7+face] = four_ninth;
    }
    else {
      coeff_out[next] += two_ninth;
      coeff_out[3]    += two_ninth;
    }
  }
  else {
    assert( face == 3);
    coeff_out[0] = -one_ninth;
    coeff_out[1] = -one_ninth;
    coeff_out[2] = -one_ninth;
    if (nodebits & 1) {
      coeff_out[4] = four_ninth;
    }
    else {
      coeff_out[0] += two_ninth;
      coeff_out[1] += two_ninth;
    }
    if (nodebits & 2) {
      coeff_out[5] = four_ninth;
    }
    else {
      coeff_out[1] += two_ninth;
      coeff_out[2] += two_ninth;
    }
    if (nodebits & 4) {
      coeff_out[6] = four_ninth;
    }
    else {
      coeff_out[2] += two_ninth;
      coeff_out[0] += two_ninth;
    }
  }
}

void TetLagrangeShape::coefficients_at_mid_elem( unsigned nodebits,
                                                 msq_std::vector<double>& coeff_out,
                                                 MsqError& err ) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  coeff_out.resize( 10 );
  coeff_out[0] = -0.125;
  coeff_out[1] = -0.125;
  coeff_out[2] = -0.125;
  coeff_out[3] = -0.125;
  if (have_node(nodebits, 4)) {
    coeff_out[4] = 0.25;
  }
  else {
    coeff_out[4] = 0.0;
    coeff_out[0] += 0.125;
    coeff_out[1] += 0.125;
  }
  if (have_node(nodebits, 5)) {
    coeff_out[5] = 0.25;
  }
  else {
    coeff_out[5] = 0.0;
    coeff_out[1] += 0.125;
    coeff_out[2] += 0.125;
  }
  if (have_node(nodebits, 6)) {
    coeff_out[6] = 0.25;
  }
  else {
    coeff_out[6] = 0.0;
    coeff_out[2] += 0.125;
    coeff_out[0] += 0.125;
  }
  if (have_node(nodebits, 7)) {
    coeff_out[7] = 0.25;
  }
  else {
    coeff_out[7] = 0.0;
    coeff_out[0] += 0.125;
    coeff_out[3] += 0.125;
  }
  if (have_node(nodebits, 8)) {
    coeff_out[8] = 0.25;
  }
  else {
    coeff_out[8] = 0.0;
    coeff_out[1] += 0.125;
    coeff_out[3] += 0.125;
  }
  if (have_node(nodebits, 9)) {
    coeff_out[9] = 0.25;
  }
  else {
    coeff_out[9] = 0.0;
    coeff_out[2] += 0.125;
    coeff_out[3] += 0.125;
  }
}

static void get_linear_derivatives( msq_std::vector<size_t>& vertices,
                                    msq_std::vector<double>& derivs )
{
  vertices.resize(4);
  derivs.resize(12);
  std::vector<size_t>::iterator v = vertices.begin();
  std::vector<double>::iterator d = derivs.begin();
  
  *v = 0; ++v;
  *d = 1.0; ++d;
  *d = 0.0; ++d;
  *d = 0.0; ++d;
  
  *v = 1; ++v;
  *d = 0.0; ++d;
  *d = 1.0; ++d;
  *d = 0.0; ++d;
  
  *v = 2; ++v;
  *d = 0.0; ++d;
  *d = 0.0; ++d;
  *d = 1.0; ++d;
  
  *v = 3; ++v;
  *d =-1.0; ++d;
  *d =-1.0; ++d;
  *d =-1.0; ++d;
}

static const unsigned edges[][2] = { { 0, 1 },
                                     { 1, 2 },
                                     { 2, 0 },
                                     { 0, 3 },
                                     { 1, 3 },
                                     { 2, 3 } };

void TetLagrangeShape::derivatives_at_corner( unsigned corner,
                                              unsigned nodebits,
                                              msq_std::vector<size_t>& vertices,
                                              msq_std::vector<double>& derivs,
                                              MsqError& err) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
    // begin with derivatives for linear tetrahedron
  get_linear_derivatives( vertices, derivs );
  
    // adjust for the presence of mid-edge nodes
  switch (corner) {
    case 0:
      if (have_node(nodebits, 4)) {
        vertices.push_back( 4 );
        derivs.push_back( 0.0 );
        derivs.push_back( 4.0 );
        derivs.push_back( 0.0 );
        derivs[1] -= 2.0;
        derivs[4] -= 2.0;
      }
      if (have_node(nodebits, 6)) {
        vertices.push_back( 6 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 4.0 );
        derivs[2] -= 2.0;
        derivs[8] -= 2.0;
      }
      if (have_node(nodebits, 7)) {
        vertices.push_back( 7 );
        derivs.push_back( -4.0 );
        derivs.push_back( -4.0 );
        derivs.push_back( -4.0 );
        derivs[ 0] += 2.0;
        derivs[ 1] += 2.0;
        derivs[ 2] += 2.0;
        derivs[ 9] += 2.0;
        derivs[10] += 2.0;
        derivs[11] += 2.0;
      }
      break;

    case 1:
      if (have_node(nodebits, 4)) {
        vertices.push_back( 4 );
        derivs.push_back( 4.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[ 0] -= 2.0;
        derivs[ 3] -= 2.0;
      }
      if (have_node(nodebits, 5)) {
        vertices.push_back(5);
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 4.0 );
        derivs[ 5] -= 2.0;
        derivs[ 8] -= 2.0;
      }
      if (have_node(nodebits, 8)) {
        vertices.push_back( 8 );
        derivs.push_back( -4.0 );
        derivs.push_back( -4.0 );
        derivs.push_back( -4.0 );
        derivs[ 3] += 2.0;
        derivs[ 4] += 2.0;
        derivs[ 5] += 2.0;
        derivs[ 9] += 2.0;
        derivs[10] += 2.0;
        derivs[11] += 2.0;
      }
      break;
  
    case 2:
      if (have_node(nodebits, 5)) {
        vertices.push_back( 5 );
        derivs.push_back( 0.0 );
        derivs.push_back( 4.0 );
        derivs.push_back( 0.0 );
        derivs[ 4] -= 2.0;
        derivs[ 7] -= 2.0;
      }
      if (have_node(nodebits, 6)) {
        vertices.push_back( 6 );
        derivs.push_back( 4.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[ 0] -= 2.0;
        derivs[ 6] -= 2.0;
      }
      if (have_node(nodebits, 9)) {
        vertices.push_back( 9 );
        derivs.push_back( -4.0 );
        derivs.push_back( -4.0 );
        derivs.push_back( -4.0 );
        derivs[ 6] += 2.0;
        derivs[ 7] += 2.0;
        derivs[ 8] += 2.0;
        derivs[ 9] += 2.0;
        derivs[10] += 2.0;
        derivs[11] += 2.0;
      }
      break;
  
    case 3:
      if (have_node(nodebits, 7)) {
        vertices.push_back( 7 );
        derivs.push_back( 4.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[ 0] -= 2.0;
        derivs[ 9] -= 2.0;
      }
      if (have_node(nodebits, 8)) {
        vertices.push_back( 8 );
        derivs.push_back( 0.0 );
        derivs.push_back( 4.0 );
        derivs.push_back( 0.0 );
        derivs[ 4] -= 2.0;
        derivs[10] -= 2.0;
      }
      
      if (have_node(nodebits, 9)) {
        vertices.push_back( 9 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 4.0 );
        derivs[ 8]-= 2.0;
        derivs[11]-= 2.0;
      }
      break;
  }
}
  
void TetLagrangeShape::derivatives_at_mid_edge( unsigned edge,
                                              unsigned nodebits,
                                              msq_std::vector<size_t>& vertices,
                                              msq_std::vector<double>& derivs,
                                              MsqError& err) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  vertices.clear();
  derivs.clear();
  
  switch (edge) {
    case 0:
      vertices.push_back( 0 );
      derivs.push_back( 1.0 );
      derivs.push_back( 0.0 );
      derivs.push_back( 0.0 );
      
      vertices.push_back( 1 );
      derivs.push_back( 0.0 );
      derivs.push_back( 1.0 );
      derivs.push_back( 0.0 );
      
   
      if (have_node(nodebits,5) && have_node(nodebits,6)) {
        vertices.push_back( 2 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back(-1.0 );
      }
      else if (!have_node(nodebits,5) && !have_node(nodebits,6)) {
        vertices.push_back( 2 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 1.0 );
      }
      
      if (!have_node(nodebits, 7) && !have_node(nodebits, 8)) {
        vertices.push_back( 3 );
        derivs.push_back( -1.0 );
        derivs.push_back( -1.0 );
        derivs.push_back( -1.0 );
      }
      else if (have_node(nodebits, 7) && have_node(nodebits, 8)) {
        vertices.push_back( 3 );
        derivs.push_back( 1.0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 1.0 );
      }
      
      if (have_node(nodebits, 4)) {
        vertices.push_back( 4 );
        derivs.push_back( 2.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs[0] -= 1.0;
        derivs[1] -= 1.0;
        derivs[3] -= 1.0;
        derivs[4] -= 1.0;
      }
      
      if (have_node(nodebits, 5)) {
        vertices.push_back( 5 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs[ 5] -= 1.0;
      }
      
      if (have_node(nodebits, 6)) {
        vertices.push_back( 6 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs[ 2] -= 1.0;
      }
      
      if (have_node(nodebits, 7)) {
        vertices.push_back( 7 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs[0] += 1.0;
        derivs[1] += 1.0;
        derivs[2] += 1.0;
      }
      
      if (have_node(nodebits, 8)) {
        vertices.push_back( 8 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs[3] += 1.0;
        derivs[4] += 1.0;
        derivs[5] += 1.0;
      }
      break;
    
    case 1:
      vertices.push_back( 1 );
      derivs.push_back( 0.0 );
      derivs.push_back( 1.0 );
      derivs.push_back( 0.0 );
      
      vertices.push_back( 2 );
      derivs.push_back( 0.0 );
      derivs.push_back( 0.0 );
      derivs.push_back( 1.0 );
   
      if (have_node(nodebits,4) && have_node(nodebits,6)) {
        vertices.push_back( 0 );
        derivs.push_back(-1.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
      }
      else if (!have_node(nodebits,4) && !have_node(nodebits,6)) {
        vertices.push_back( 0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
      }
      
      if (!have_node(nodebits, 8) && !have_node(nodebits, 9)) {
        vertices.push_back( 3 );
        derivs.push_back( -1.0 );
        derivs.push_back( -1.0 );
        derivs.push_back( -1.0 );
      }
      else if (have_node(nodebits, 8) && have_node(nodebits, 9)) {
        vertices.push_back( 3 );
        derivs.push_back( 1.0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 1.0 );
      }
      
      if (have_node(nodebits, 4)) {
        vertices.push_back( 4 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[0] -= 1.0;
      }
      
      if (have_node(nodebits, 5)) {
        vertices.push_back( 5 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 2.0 );
        derivs[1] -= 1.0;
        derivs[2] -= 1.0;
        derivs[4] -= 1.0;
        derivs[5] -= 1.0;
      }
      
      if (have_node(nodebits, 6)) {
        vertices.push_back( 6 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[3] -= 1.0;
      }
      
      if (have_node(nodebits, 8)) {
        vertices.push_back( 8 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs[0] += 1.0;
        derivs[1] += 1.0;
        derivs[2] += 1.0;
      }
      
      if (have_node(nodebits, 9)) {
        vertices.push_back( 9 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs[3] += 1.0;
        derivs[4] += 1.0;
        derivs[5] += 1.0;
      }
      break;
      
    case 2:
      vertices.push_back( 0 );
      derivs.push_back( 1.0 );
      derivs.push_back( 0.0 );
      derivs.push_back( 0.0 );
      
      vertices.push_back( 2 );
      derivs.push_back( 0.0 );
      derivs.push_back( 0.0 );
      derivs.push_back( 1.0 );
   
      if (have_node(nodebits,4) && have_node(nodebits,5)) {
        vertices.push_back( 1 );
        derivs.push_back( 0.0 );
        derivs.push_back(-1.0 );
        derivs.push_back( 0.0 );
      }
      if (!have_node(nodebits,4) && !have_node(nodebits,5)) {
        vertices.push_back( 1 );
        derivs.push_back( 0.0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 0.0 );
      }
      
      if (!have_node(nodebits, 7) && !have_node(nodebits, 9)) {
        vertices.push_back( 3 );
        derivs.push_back( -1.0 );
        derivs.push_back( -1.0 );
        derivs.push_back( -1.0 );
      }
      else if (have_node(nodebits, 7) && have_node(nodebits, 9)) {
        vertices.push_back( 3 );
        derivs.push_back( 1.0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 1.0 );
      }
      
      if (have_node(nodebits, 4)) {
        vertices.push_back( 4 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs[ 1] -= 1.0;
      }
      
      if (have_node(nodebits, 5)) {
        vertices.push_back( 5 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs[ 4] -= 1.0;
      }
      
      if (have_node(nodebits, 6)) {
        vertices.push_back( 6 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs[0] -= 1.0;
        derivs[2] -= 1.0;
        derivs[3] -= 1.0;
        derivs[5] -= 1.0;
      }
       
      if (have_node(nodebits, 7)) {
        vertices.push_back( 7 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs[0] += 1.0;
        derivs[1] += 1.0;
        derivs[2] += 1.0;
      }
      
      if (have_node(nodebits, 9)) {
        vertices.push_back( 9 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs[3] += 1.0;
        derivs[4] += 1.0;
        derivs[5] += 1.0;
      }
      break;
    
    case 3:
      vertices.push_back( 0 );
      derivs.push_back( 1.0 );
      derivs.push_back( 0.0 );
      derivs.push_back( 0.0 );
      
      vertices.push_back( 3 );
      derivs.push_back( -1.0 );
      derivs.push_back( -1.0 );
      derivs.push_back( -1.0 );
      
      if (!have_node( nodebits, 4 ) && !have_node( nodebits, 8 )) {
        vertices.push_back( 1 );
        derivs.push_back( 0.0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 0.0 );
      }
      else if (have_node( nodebits, 4 ) && have_node( nodebits, 8 )) {
        vertices.push_back( 1 );
        derivs.push_back( 0.0 );
        derivs.push_back(-1.0 );
        derivs.push_back( 0.0 );
      }
      
      if (!have_node( nodebits, 6 ) && !have_node( nodebits, 9 )) {
        vertices.push_back( 2 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 1.0 );
      }
      if (have_node( nodebits, 6 ) && have_node( nodebits, 9 )) {
        vertices.push_back( 2 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back(-1.0 );
      }

      if (have_node( nodebits, 4 )) {
        vertices.push_back( 4 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs[1] -= 1.0;
      }
      
      if (have_node( nodebits, 6 )) {
        vertices.push_back( 6 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs[2] -= 1.0;
      }
      
      if (have_node( nodebits, 7 )) {
        vertices.push_back( 7 );
        derivs.push_back(  0.0 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs[1] += 1.0;
        derivs[2] += 1.0;
        derivs[4] += 1.0;
        derivs[5] += 1.0;
      }
      
      if (have_node( nodebits, 8 )) {
        vertices.push_back( 8 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs[4] -= 1.0;
      }
      if (have_node( nodebits, 9 )) {
        vertices.push_back( 9 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs[5] -= 1.0;
      }
      break;
    
    case 4:
      vertices.push_back( 1 );
      derivs.push_back( 0.0 );
      derivs.push_back( 1.0 );
      derivs.push_back( 0.0 );
      
      vertices.push_back( 3 );
      derivs.push_back( -1.0 );
      derivs.push_back( -1.0 );
      derivs.push_back( -1.0 );

      if (!have_node( nodebits, 4 ) && !have_node( nodebits, 7 )) {
        vertices.push_back( 0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
      }
      else if (have_node( nodebits, 4 ) && have_node( nodebits, 7 )) {
        vertices.push_back( 0 );
        derivs.push_back(-1.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
      }
      
      if (!have_node( nodebits, 5 ) && !have_node( nodebits, 9 )) {
        vertices.push_back( 2 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 1.0 );
      }
      else if (have_node( nodebits, 5 ) && have_node( nodebits, 9 )) {
        vertices.push_back( 2 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back(-1.0 );
      }

      if (have_node( nodebits, 4 )) {
        vertices.push_back( 4 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[0] -= 1.0;
      }
      
      if (have_node( nodebits, 5 )) {
        vertices.push_back( 5 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs[2] -= 1.0;
      }
      
      if (have_node( nodebits, 8 )) {
        vertices.push_back( 8 );
        derivs.push_back( -2.0 );
        derivs.push_back(  0.0 );
        derivs.push_back( -2.0 );
        derivs[0] += 1.0;
        derivs[2] += 1.0;
        derivs[3] += 1.0;
        derivs[5] += 1.0;
      }
      
      if (have_node( nodebits, 7 )) {
        vertices.push_back( 7 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[3] -= 1.0;
      }
      if (have_node( nodebits, 9 )) {
        vertices.push_back( 9 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs[5] -= 1.0;
      }
      break;
    
    case 5:
      vertices.push_back( 2 );
      derivs.push_back( 0.0 );
      derivs.push_back( 0.0 );
      derivs.push_back( 1.0 );
      
      vertices.push_back( 3 );
      derivs.push_back( -1.0 );
      derivs.push_back( -1.0 );
      derivs.push_back( -1.0 );

      if (!have_node( nodebits, 6 ) && !have_node( nodebits, 7 )) {
        vertices.push_back( 0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
      }
      else if (have_node( nodebits, 6 ) && have_node( nodebits, 7 )) {
        vertices.push_back( 0 );
        derivs.push_back(-1.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
      }
      
      if (!have_node( nodebits, 5 ) && !have_node( nodebits, 8 )) {
        vertices.push_back( 1 );
        derivs.push_back( 0.0 );
        derivs.push_back( 1.0 );
        derivs.push_back( 0.0 );
      }
      else if (have_node( nodebits, 5 ) && have_node( nodebits, 8 )) {
        vertices.push_back( 1 );
        derivs.push_back( 0.0 );
        derivs.push_back(-1.0 );
        derivs.push_back( 0.0 );
      }
      
      if (have_node( nodebits, 5 )) {
        vertices.push_back( 5 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs[1] -= 1.0;
      }

      if (have_node( nodebits, 6 )) {
        vertices.push_back( 6 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[0] -= 1.0;
      }
      
      if (have_node( nodebits, 9 )) {
        vertices.push_back( 9 );
        derivs.push_back( -2.0 );
        derivs.push_back( -2.0 );
        derivs.push_back(  0.0 );
        derivs[0] += 1.0;
        derivs[1] += 1.0;
        derivs[3] += 1.0;
        derivs[4] += 1.0;
      }
      
      if (have_node( nodebits, 7 )) {
        vertices.push_back( 7 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs.push_back( 0.0 );
        derivs[3] -= 1.0;
      }
      if (have_node( nodebits, 8 )) {
        vertices.push_back( 8 );
        derivs.push_back( 0.0 );
        derivs.push_back( 2.0 );
        derivs.push_back( 0.0 );
        derivs[4] -= 1.0;
      }
      break;
  }
}

// Derivatives of coefficients for higher-order nodes

const double ft = 4.0/3.0;

const double ho_dr[6][4] = { { ft, ft, 0., ft },
                             { 0., 0., 0., 0. },
                             { 0., ft, ft, ft },
                             { 0., ft, 0.,-ft },
                             {-ft,-ft, 0.,-ft },
                             { 0.,-ft,-ft,-ft } };

const double ho_ds[6][4] = { { ft, 0., ft, ft },
                             { 0., ft, ft, ft },
                             { 0., 0., 0., 0. },
                             {-ft, 0.,-ft,-ft },
                             { 0., 0., ft,-ft },
                             { 0.,-ft,-ft,-ft } };

const double ho_dt[6][4] = { { 0., 0., 0., 0. },
                             { ft, ft, 0., ft },
                             { ft, 0., ft, ft },
                             {-ft, 0.,-ft,-ft },
                             {-ft,-ft, 0.,-ft },
                             { ft, 0., 0.,-ft } };


void TetLagrangeShape::derivatives_at_mid_face( unsigned face,
                                              unsigned nodebits,
                                              msq_std::vector<size_t>& vertices,
                                              msq_std::vector<double>& derivs,
                                              MsqError& err) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
    // begin with derivatives for linear tetrahedron
  get_linear_derivatives( vertices, derivs );
  
  for (unsigned i = 0; i < 6; ++i) 
    if (nodebits & (1<<i)) {
      vertices.push_back( i+4 );
      derivs.push_back( ho_dr[i][face] );
      derivs.push_back( ho_ds[i][face] );
      derivs.push_back( ho_dt[i][face] );
      int j = 3*edges[i][0];
      derivs[j  ] -= 0.5*ho_dr[i][face];
      derivs[j+1] -= 0.5*ho_ds[i][face];
      derivs[j+2] -= 0.5*ho_dt[i][face];
      j = 3*edges[i][1];
      derivs[j  ] -= 0.5*ho_dr[i][face];
      derivs[j+1] -= 0.5*ho_ds[i][face];
      derivs[j+2] -= 0.5*ho_dt[i][face];
    }
}

void TetLagrangeShape::derivatives_at_mid_elem( 
                                              unsigned nodebits,
                                              msq_std::vector<size_t>& vertices,
                                              msq_std::vector<double>& derivs,
                                              MsqError& err) const
{
  if (nodebits >= (1u << 6)) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
                                    
  bool corners[4] = { false, false, false, false };
  double corner_vals[4][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
  
  vertices.clear();
  derivs.clear();
  
  for (unsigned i = 4;  i < 10; ++i) {
    int sign, zero;
    if (i < 7) {
      sign = 1;
      zero = (i - 2) % 3;
    }
    else {
      sign = -1;
      zero = (i - 7);
    }
    
    if (have_node( nodebits, i )) {
      vertices.push_back( i );
      int n = derivs.size();
      derivs.resize( n+3, (double)sign );
      derivs[n+zero] = 0.0;
    }
    else {
      for (unsigned j = 0; j < 2; ++j) {
        int corner = edges[i-4][j];
        int v1 = (zero + 1) % 3;
        int v2 = (zero + 2) % 3;
        corners[corner] = true;
        corner_vals[corner][v1] += 0.5*sign;
        corner_vals[corner][v2] += 0.5*sign;
      }
    }
  }
  
  for (unsigned i = 0; i < 4; ++i)
    if (corners[i]) {
      vertices.push_back(i);
      derivs.push_back(corner_vals[i][0]);
      derivs.push_back(corner_vals[i][1]);
      derivs.push_back(corner_vals[i][2]);
    }
}
    


} // namespace Mesquite
