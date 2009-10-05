/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_INSTRUCTION_HPP
#define MSQ_INSTRUCTION_HPP

#include "Mesquite.hpp"

#include <string>

namespace MESQUITE_NS {

class Mesh;
class ParallelMesh;
class MeshDomain;
class MsqError;
class Settings;

class MESQUITE_EXPORT Instruction
{
  public:
  
    virtual ~Instruction() {}
  
    virtual double loop_over_mesh( Mesh* mesh, 
                                   MeshDomain* domain, 
                                   const Settings* settings,
                                   MsqError& err ) = 0;

    virtual double loop_over_mesh( ParallelMesh* mesh, 
                                   MeshDomain* domain, 
                                   const Settings* settings,
                                   MsqError& err );

    virtual std::string get_name() const = 0;
};

} // namespace Mesquite

#endif
