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


/** \file TargetWriter.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_WRITER_HPP
#define MSQ_TARGET_WRITER_HPP

#include "Mesquite.hpp"
#include "Instruction.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS {

class TargetCalculator;
class WeightCalculator;

/** Save target matrices in tag data 
 *
 * Store element target matrices in tag data.  The stored target
 * matrices may be retreived using the TargetReader class.
 */
class TargetWriter : public Instruction
{
public:
  TargetWriter( TargetCalculator* tc,
                WeightCalculator* wc = 0,
                msq_std::string target_base_name = "MSQ_TARGET_MATRIX",
                msq_std::string weight_base_name = "MSQ_TARGET_WEIGHT") ;
                
  virtual ~TargetWriter();
  
  double loop_over_mesh( Mesh* mesh, MeshDomain*, const Settings*, MsqError& );

  msq_std::string get_name() const;

private:

  TagHandle get_target_tag( unsigned dimension,  unsigned count, 
                            Mesh* mesh, MsqError& err );
  TagHandle get_weight_tag( unsigned count, Mesh* mesh, MsqError& err );
  TagHandle get_tag_handle( const msq_std::string& base_name,
                            unsigned num_dbl, Mesh* mesh, MsqError& err );

  TargetCalculator* targetCalc;
  WeightCalculator* weightCalc;
  
  msq_std::string targetName, weightName;
  msq_std::vector<TagHandle> targetTags, weightTags;
};


} // namespace Mesquite

#endif
