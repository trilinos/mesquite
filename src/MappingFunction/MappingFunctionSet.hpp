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


/** \file MappingFunctionSet.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MAPPING_FUNCTION_SET_HPP
#define MSQ_MAPPING_FUNCTION_SET_HPP

#include "Mesquite.hpp"

namespace Mesquite {

class MappingFunction;
class MappingFunction2D;
class MappingFunction3D;

/**\brief Class describing mapping functions for each element type
 *
 * An interface for a class that manages the list of mapping functions
 *for each element type.
 */
class MappingFunctionSet
{
public:
  virtual ~MappingFunctionSet() {}
    /** Get mapping function for specified element type */
  virtual const MappingFunction* get_function( EntityTopology elem_type ) const = 0;
    /** Get mapping function for surface element */
  virtual const MappingFunction2D* get_surf_function( EntityTopology elem_type ) const = 0;
    /** Get mapping function for volume element */
  virtual const MappingFunction3D* get_vol_function( EntityTopology elem_type ) const = 0;
};


} // namespace Mesquite

#endif
