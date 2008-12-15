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


/** \file LinearFunctionSet.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LinearFunctionSet.hpp"

#include "LinearTriangle.hpp"
#include "LinearQuadrilateral.hpp"
#include "LinearTetrahedron.hpp"
#include "LinearPyramid.hpp"
#include "LinearPrism.hpp"
#include "LinearHexahedron.hpp"

namespace Mesquite {

static const LinearTriangle tri;
static const LinearQuadrilateral quad;
static const LinearTetrahedron tet;
static const LinearPyramid pyr;
static const LinearPrism prism;
static const LinearHexahedron hex;

static const MappingFunction* func_array[MIXED] =
 { 0, 0, 0, 0, 0, 0, 0, 0, &tri, &quad, 0, &tet, &hex, &prism, &pyr, 0 };

static const MappingFunction2D* func_array_2d[MIXED] =
 { 0, 0, 0, 0, 0, 0, 0, 0, &tri, &quad, 0, 0, 0, 0, 0, 0 };

static const MappingFunction3D* func_array_3d[MIXED] =
 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &tet, &hex, &prism, &pyr, 0 };

const MappingFunction* LinearFunctionSet::get_function( EntityTopology type ) const
  { return func_array[type]; }

const MappingFunction2D* LinearFunctionSet::get_surf_function( EntityTopology type ) const
  { return func_array_2d[type]; }

const MappingFunction3D* LinearFunctionSet::get_vol_function( EntityTopology type ) const
  { return func_array_3d[type]; }

} // namespace Mesquite
