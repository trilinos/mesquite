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


/** \file SamplePoints.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SAMPLE_POINTS_HPP
#define MSQ_SAMPLE_POINTS_HPP

#include "Mesquite.hpp"
#include <assert.h>

namespace Mesquite {

/**\brief Topological description of element sample points.
 *
 * This class provides a set of static methods for working with
 * topologically defined sample points and a method for storing
 * per-element-type sample point specifications.
 *
 * Sample points are defined by the topological dimension of the 
 * entity they lie at the center of.  For example:
 * - 0 corner
 * - 1 mid-edge
 * - 2 mid-face
 *
 * The total sample points for a given element type is stored
 * in a single numberical value, one bit per topological sub-entity
 * dimension.  The position of each bit is the dimension of the 
 * sub-entity, where the least signficant bit is zero.  This representation
 * is used primarily for internal member data for this class.  It is
 * recommended that other code treat such values as opaque and use
 * dimension-based methods whenever possible.
 */
struct SampleLoc { char dim; char num; };
typedef const SampleLoc* const (*loc_table_ptr)[16];
class SamplePoints
{
public:
  
  /**\brief Get number of sample points for an element type
   *
   * Given a bit mask describing the topological locations of
   * sample points in an element, return the total count of
   * sample points.
   */
  static unsigned num_sample_points( EntityTopology element_type, 
                                     unsigned sample_topology_bits );
  
  /**\brief Get the sample point number given the topological location.
   *
   * Get the sample point number for an element type.
   *
   *\param element_type  The element topology
   *\param sample_topology_bits Bit set describing all element sample points
   *\param entity_dimension The dimension of the entity for which
   *                           the sample point is the topological midpoint
   *                           (e.g. 1->mid_edge)
   *\param entity_number The number of the entity of dimension entity_dimensino
   *                     in the canonical order for the element topology.
   *\return  A number in the range [0,num_sample_points], or (unsigned)-1
   *         if the there are no sample points for the specified sub-entity
   *         dimension.
   */
  static unsigned sample_number_from_location( EntityTopology element_type,
                                           unsigned sample_topology_bits,
                                           unsigned entity_dimension,
                                           unsigned entity_number );
 
  /**\brief Initialize sample points for each element type to either none 
    * or the default for that element type. */
  SamplePoints( bool do_default = false );
 
  /**\brief Initialize sample points for each element type. */
  SamplePoints( bool corners,
                bool mid_edge,
                bool mid_face,
                bool mid_elem);
                
  /**\brief Clear entries for all ElementTopologies. */
  void clear_all();

  /**\brief Get number of sample points for an element type
   *
   * Given a bit mask describing the topological locations of
   * sample points in an element, return the total count of
   * sample points.
   */
  inline unsigned num_sample_points( EntityTopology type ) const;
  
  /**\brief Test if element will be sampled at the specified topological location
   *
   * Test if the specified element topology will be sampled at the
   * mid-point of sub-entities with the specified dimension.
   */
  inline bool will_sample_at( EntityTopology type, unsigned dimension ) const;
  
  /**\brief Sample elements at the specified topological location
   *
   * Sample the specified element topology at the
   * mid-point of sub-entities with the specified dimension.
   */
  inline void sample_at( EntityTopology type, unsigned dimension );
  
  /**\brief Do not sample elements at the specified topological location
   *
   * Do not sample the specified element topology at the
   * mid-point of sub-entities with the specified dimension.
   */
  inline void dont_sample_at( EntityTopology type, unsigned dimension );
  
  /**\brief Get bit set describing element sample points */
  inline unsigned get_sample_topology_bits( EntityTopology type ) const;
 
  /**\brief Get the sample point number given the topological location.
   *
   * Get the sample point number for an element type.
   *
   *\param element_type  The element topology
   *\param entity_dimension The dimension of the entity for which
   *                           the sample point is the topological midpoint
   *                           (e.g. 1->mid_edge)
   *\param entity_number The number of the entity of dimension entity_dimensino
   *                     in the canonical order for the element topology.
   *\return  A number in the range [0,num_sample_points], or (unsigned)-1
   *         if the there are no sample points for the specified sub-entity
   *         dimension.
   */
  inline unsigned sample_number_from_location( EntityTopology element_type,
                                           unsigned entity_dimension,
                                           unsigned entity_number ) const;
  
  /**\brief Set bit set describing element sample points */
  inline void set_sample_topology_bits( EntityTopology type, unsigned bits );
  
private:

  unsigned char sampleBits[Mesquite::MIXED];
  static loc_table_ptr locationTable;
};
  
  
unsigned SamplePoints::get_sample_topology_bits( EntityTopology type ) const
  { return sampleBits[type]; }

void SamplePoints::set_sample_topology_bits( EntityTopology type, unsigned bits )
  { sampleBits[type] = bits; }
 
unsigned SamplePoints::num_sample_points( EntityTopology type ) const
  { return num_sample_points( type, get_sample_topology_bits( type ) ); }
  
bool SamplePoints::will_sample_at( EntityTopology type, unsigned dimension ) const
  { return 0 != (sampleBits[type] & (1 << dimension)); }
  
void SamplePoints::sample_at( EntityTopology type, unsigned dimension )
  { sampleBits[type] |= (1 << dimension ); }
  
void SamplePoints::dont_sample_at( EntityTopology type, unsigned dimension )
  { sampleBits[type] &= ~(1 << dimension ); }
  
unsigned SamplePoints::sample_number_from_location( EntityTopology element_type,
                                                    unsigned entity_dimension,
                                                    unsigned entity_number ) const
{
  return sample_number_from_location( element_type,
                                      get_sample_topology_bits( element_type ),
                                      entity_dimension, entity_number );
}

} // namespace Mesquite

#endif
