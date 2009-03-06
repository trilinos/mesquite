/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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


/** \file NodeSet.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_NODE_SET_HPP
#define MSQ_NODE_SET_HPP

#include "Mesquite.hpp"
#include <assert.h>
#include <iosfwd>
#include "Bits.hpp"
#include "Sample.hpp"

namespace Mesquite {

/** Utility class for storing one boolean mark/flag for each node in an element */
class NodeSet {
  private:
    typedef unsigned BitSet;
    BitSet bits; //!< The data, one bit for each possible node location
    
    //! Return a bool value indicating the state of the specified bit.
    inline bool bit_to_bool( unsigned num ) const
      { return static_cast<bool>( (bits >> num) & 1u ); }
    
    //! Set the specified bit to 1
    inline void set_bit( unsigned num )
      { bits |= (1u << num); }
    
    //! Clear the specified bit (set it to zero)
    inline void clear_bit( unsigned num )
      { bits &= ~(1u << num); }
    
    //! Count the number of non-zero bits with less significance than
    //! the specified bit.  (Mask out all bits except those of before
    //! the specified position and then count the remaining bits.)
    unsigned num_before_bit( unsigned position ) const
      { return popcount(bits & ~((~0u) << position)); }
  
  public:
    
    //! Misc constants.  Only NUM_CORNER_BITS, NUM_EDGE_BITS, and
    //! NUM_REGION_BITS should be modified.  All other contants are
    //! a function of those three values and the size of the bit storage.
    enum { 
      NUM_TOTAL_BITS = 8*sizeof(BitSet),
      MSB_POS = NUM_TOTAL_BITS - 1,
      //! Maximum number of corner nodes.
      NUM_CORNER_BITS = 8, 
      //! Maximum number of mid-edge nodes
      NUM_EDGE_BITS = 16,  
      //! Maximum number of mid-volume nodes
      NUM_REGION_BITS = 1, 
      //! Maximum number of mid-face nodes
      NUM_FACE_BITS = NUM_TOTAL_BITS - (NUM_CORNER_BITS + NUM_EDGE_BITS + NUM_REGION_BITS),   
      
      //! LSB of corner node storage
      CORNER_OFFSET = 0,
      //! LSB of mid-edgee storage
      EDGE_OFFSET = CORNER_OFFSET + NUM_CORNER_BITS,
      //! LSB of mid-face storage
      FACE_OFFSET =   EDGE_OFFSET + NUM_EDGE_BITS,
      //! LSB of mid-region storage
      REGION_OFFSET = FACE_OFFSET + NUM_FACE_BITS,
    
      //! Bit mask for all non-corner bits
      MID_NODE_MASK = (~0u) << NUM_CORNER_BITS,
      //! Bit mask for all corner bits
      CORNER_MASK = ~MID_NODE_MASK,
      //! Mid-region mask
      REGION_MASK = ~((~0u) >> NUM_REGION_BITS),
      //! Bit mask for all mid-edge nodes
      EDGE_MASK = (MID_NODE_MASK << (NUM_TOTAL_BITS - EDGE_OFFSET - NUM_EDGE_BITS)) >> (NUM_TOTAL_BITS - EDGE_OFFSET - NUM_EDGE_BITS),
      //! Bit mask for all mid-face nodes
      FACE_MASK = ~(CORNER_MASK|EDGE_MASK|REGION_MASK),
    };
    
    NodeSet() : bits(0u) {}
    
    //! Clear all values
    void clear() { bits = 0u; }
    
    //! Set all marks/flags
    void set_all_nodes() { bits = ~0u; }
    
    //! Number of marked/flagged nodes
    unsigned num_nodes() const
      { return popcount(bits); }
      
    //! Check if any mid-nodes (higher-order nodes) are flaged
    BitSet have_any_mid_node() const
      { return bits & MID_NODE_MASK; }
    //! Check if any corner nodes are present
    BitSet have_any_corner_node() const
      { return bits & CORNER_MASK; }
    //! Check if any mid-edge nodes are present
    BitSet have_any_mid_edge_node() const
      { return bits & EDGE_MASK; }
    //! Check if any mid-face nodes are present
    BitSet have_any_mid_face_node() const
      { return bits & FACE_MASK; }
    //! Check if any mid-region nodes are present
    BitSet have_any_mid_region_node() const
      { return bits & REGION_MASK; }
      
    //! Set all mid-nodes
    void set_all_mid_nodes()
      { bits |= MID_NODE_MASK; }
    //! Set all corner nodes
    void set_all_corner_nodes()
      { bits |= CORNER_MASK; }
    //! Set all mid-edge nodes
    void set_all_mid_edge_nodes()
      { bits |= EDGE_MASK; }
    //! Set all mid-face nodes
    void set_all_mid_face_nodes()
      { bits |= FACE_MASK; }
    //! Set all mid-region nodes
    void set_all_mid_region_nodes()
      { bits |= REGION_MASK; }
      
    //! Clear all mid-nodes
    void clear_all_mid_nodes()
      { bits &= ~MID_NODE_MASK; }
    //! Clear all corner nodes
    void clear_all_corner_nodes()
      { bits &= ~CORNER_MASK; }
    //! Clear all mid-edge nodes
    void clear_all_mid_edge_nodes()
      { bits &= ~EDGE_MASK; }
    //! Clear all mid-face nodes
    void clear_all_mid_face_nodes()
      { bits &= ~FACE_MASK; }
    //! Clear all mid-region nodes
    void clear_all_mid_region_nodes()
      { bits &= ~REGION_MASK; }
    
    //! Position of a corner node in the list
    static unsigned corner_node_position( unsigned num ) 
      { assert(num < NUM_CORNER_BITS); return num + CORNER_OFFSET; }
    //! Position of a mid-edge node in the list
    static unsigned mid_edge_node_position( unsigned num ) 
      { assert(num < NUM_EDGE_BITS); return num + EDGE_OFFSET; }
    //! Position of a mid-face node in the list
    static unsigned mid_face_node_position( unsigned num ) 
      { assert(num < NUM_FACE_BITS); return num + FACE_OFFSET; }
    //! Position of a mid-region node in the list
    static unsigned mid_region_node_position( unsigned num = 0 ) 
      { assert(num < NUM_REGION_BITS); return num + REGION_OFFSET; }
    //! Position of a node in the list
    static unsigned position( Sample sample ) {
      switch (sample.dimension) { 
        case 0: return corner_node_position(sample.number); 
        case 1: return mid_edge_node_position(sample.number);
        case 2: return mid_face_node_position(sample.number);
        case 3: return mid_region_node_position(sample.number);
#ifndef NDEBUG
        default: assert(0); return ~0u;
#endif
      }
    }
    
    //! Mark/flag corner node
    void set_corner_node( unsigned num ) 
      { set_bit( corner_node_position(num) ); }
    //! Mark/flag mid-edge node
    void set_mid_edge_node( unsigned num ) 
      { set_bit( mid_edge_node_position(num) ); }
    //! Mark/flag mid-face node
    void set_mid_face_node( unsigned num ) 
      { set_bit( mid_face_node_position(num) ); }
    //! Mark/flag mid-region node
    void set_mid_region_node( unsigned num = 0 ) 
      { set_bit( mid_region_node_position(num) ); }
    //! Mark/flag node
    void set_node( Sample sample ) {
      switch (sample.dimension) { 
        case 0: set_corner_node(sample.number); break;
        case 1: set_mid_edge_node(sample.number); break;
        case 2: set_mid_face_node(sample.number); break;
        case 3: set_mid_region_node(sample.number); break;
        default: assert(0);
      }
    }
    
    //! un-mark (clear flag for) corner node
    void clear_corner_node( unsigned num ) 
      { clear_bit( corner_node_position(num) ); }
    //! un-mark (clear flag for) mid-edge node
    void clear_mid_edge_node( unsigned num ) 
      { clear_bit( mid_edge_node_position(num) ); }
    //! un-mark (clear flag for) mid-face node
    void clear_mid_face_node( unsigned num ) 
      { clear_bit( mid_face_node_position(num) ); }
    //! un-mark (clear flag for) mid-region node
    void clear_mid_region_node( unsigned num = 0 ) 
      { clear_bit( mid_region_node_position(num) ); }
    //! un-mark (clear flag for) node
    void clear_node( Sample sample ) {
      switch (sample.dimension) { 
        case 0: clear_corner_node(sample.number); break;
        case 1: clear_mid_edge_node(sample.number); break;
        case 2: clear_mid_face_node(sample.number); break;
        case 3: clear_mid_region_node(sample.number); break;
        default: assert(0);
      }
    }
    
    //! Get mark/flag for corner node
    bool corner_node( unsigned num ) const
      { return bit_to_bool( corner_node_position(num) ); }
    //! Get mark/flag for mid-edge node
    bool mid_edge_node( unsigned num ) const
      { return bit_to_bool( mid_edge_node_position(num) ); }
    //! Test if two mid-edge nodes are both present
    bool both_edge_nodes( unsigned num1, unsigned num2 ) const
    { 
      BitSet b = (1<<mid_edge_node_position(num1))
               | (1<<mid_edge_node_position(num2));
      return (bits&b) == b; 
    }
    //! Get mark/flag for mid-face node
    bool mid_face_node( unsigned num ) const
      { return bit_to_bool( mid_face_node_position(num) ); }
    //! Get mark/flag for mid-region node
    bool mid_region_node( unsigned num = 0 ) const
      {  return bit_to_bool( mid_region_node_position(num) ); }
    //! Get mark/flag for node
    bool node( Sample sample ) const {
      switch (sample.dimension) { 
        case 0: return corner_node(sample.number); break;
        case 1: return mid_edge_node(sample.number); break;
        case 2: return mid_face_node(sample.number); break;
        case 3: return mid_region_node(sample.number); break;
#ifndef NDEBUG
        default: assert(0); return false;
#endif
      }
    }

    //! Get number of marked/flagged nodes preceeding a specified corner node in the list
    unsigned num_before_corner( unsigned num ) const 
      { return num_before_bit( corner_node_position(num) ); } 
    //! Get number of marked/flagged nodes preceeding a specified mid-edge node in the list
    unsigned num_before_mid_edge( unsigned num ) const 
      { return num_before_bit( mid_edge_node_position(num) ); } 
    //! Get number of marked/flagged nodes preceeding a specified mid-face node in the list
    unsigned num_before_mid_face( unsigned num ) const 
      { return num_before_bit( mid_face_node_position(num) ); } 
    //! Get number of marked/flagged nodes preceeding a specified mid-region node in the list
    unsigned num_before_mid_region( unsigned num ) const 
      { return num_before_bit( mid_region_node_position(num) ); } 
    //! Get number of marked/flagged nodes preceeding a specified node in the list
    unsigned num_before( Sample sample ) const {
      switch (sample.dimension) { 
        case 0: return num_before_corner(sample.number); break;
        case 1: return num_before_mid_edge(sample.number); break;
        case 2: return num_before_mid_face(sample.number); break;
        case 3: return num_before_mid_region(sample.number); break;
#ifndef NDEBUG
        default: assert(0); return ~0u;
#endif
      }
    }
    
    unsigned get_bits() const { return bits; }
};

//! Print bits in reverse order (least-signficant to most-significant,
//! or corner 0 to mid-region).
msq_stdio::ostream& operator<<( msq_stdio::ostream& s, NodeSet set );

} // namespace Mesquite

#endif
