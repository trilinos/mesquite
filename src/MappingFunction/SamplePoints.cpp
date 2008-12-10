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


/** \file SamplePoints.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "SamplePoints.hpp"
#include "TopologyInfo.hpp"
#include <assert.h>
#include <stdlib.h>

namespace Mesquite {

const unsigned char CORNERS = 1;
const unsigned char MIDEDGE = (CORNERS << 1);
const unsigned char MIDFACE = (CORNERS << 2);
const unsigned char MIDELEM = (CORNERS << 3);

static loc_table_ptr get_location_table();
loc_table_ptr SamplePoints::locationTable = get_location_table();

SamplePoints::SamplePoints( bool do_default )
{
  clear_all();
  if (do_default) {
    sampleBits[POLYGON      ] = CORNERS;
    sampleBits[TRIANGLE     ] = MIDFACE;
    sampleBits[QUADRILATERAL] = CORNERS;
    sampleBits[POLYHEDRON   ] = 0;
    sampleBits[TETRAHEDRON  ] = MIDELEM;
    sampleBits[HEXAHEDRON   ] = CORNERS;
    sampleBits[PRISM        ] = CORNERS;
    sampleBits[PYRAMID      ] = CORNERS;
    sampleBits[SEPTAHEDRON  ] = CORNERS;
  }
}

void SamplePoints::clear_all()
  { for (size_t i = 0; i < MIXED; sampleBits[i++] = 0); }

SamplePoints::SamplePoints( bool corners,
                            bool midedge,
                            bool midface,
                            bool midelem)
{
  size_t i;
  unsigned char bits = 0;
  if (corners) bits |= CORNERS;
  if (midedge) bits |= MIDEDGE;
  if (midface) bits |= MIDFACE;
  for (i = 0; i < TETRAHEDRON; ++i)
    sampleBits[i] = bits;
  if (midelem) bits |= MIDELEM;
  for (i = TETRAHEDRON; i < MIXED; ++i)
    sampleBits[i] = bits;
}
  
unsigned SamplePoints::num_sample_points( EntityTopology type, unsigned bits )
{
  const unsigned dim = TopologyInfo::dimension( type );
  unsigned count = !!(bits & (1 << dim));
  for (unsigned i = 0; i < dim; ++i)
    if (bits & (1 << i))
      count += TopologyInfo::adjacent( type, i );
    // special case for pyramids: we never want to evaluate at
    // the apex of the pyramid
  if (type == PYRAMID && bits&1)
    --count;
  return count;
}

/* Replace general-purpose code with inlined table lookup.
   Elimintes the 3% of CPU time spent here.
   
void SamplePoints::location_from_sample_number( EntityTopology element_type,
                                           unsigned sample_topology_bits,
                                           unsigned sample_point_number,
                                           unsigned& entity_dimension,
                                           unsigned& entity_number )
{
  const unsigned dim = TopologyInfo::dimension( element_type );
  for (unsigned i = 0; i < dim; ++i) {
    if ( sample_topology_bits & (1 << i)) {
      unsigned n =  TopologyInfo::adjacent( element_type, i );
        // Skip apex for pyramid elements
      if (element_type == PYRAMID && i == 0)
        --n;
      if (n > sample_point_number) {
        entity_number = sample_point_number;
        entity_dimension = i;
        return;
      }
      sample_point_number -= n;
    }
  }
  
  assert( sample_point_number == 0 );
  entity_dimension = dim;
  entity_number = 0;
}
*/

unsigned SamplePoints::sample_number_from_location( EntityTopology element_type,
                                           unsigned sample_topology_bits,
                                           unsigned entity_dimension,
                                           unsigned entity_number )
{
  const unsigned dim = TopologyInfo::dimension( element_type );
  assert( entity_dimension <= dim );
  if (!(sample_topology_bits & (1 << entity_dimension)))
    return (unsigned)-1; 
  
  unsigned result = (entity_dimension == dim) ? 0 : entity_number;
  for (unsigned i = 0; i < entity_dimension; ++i) 
    if (sample_topology_bits & (1 << i)) 
      result += TopologyInfo::adjacent( element_type, i );
  
    // skip apex for pyramids
  if (element_type == PYRAMID && (sample_topology_bits & 1) && result > 3)
    --result;
    
  return result;
}

/* Portions of look-up table for location_from_sample_number */
static const SampleLoc full_tri[7] = { { 0, 0 },
                                       { 0, 1 },
                                       { 0, 2 },
                                       { 1, 0 },
                                       { 1, 1 },
                                       { 1, 2 },
                                       { 2, 0 } };
static const SampleLoc tri_no_edge[4] = { { 0, 0 },
                                          { 0, 1 },
                                          { 0, 2 },
                                          { 2, 0 } };
static const SampleLoc full_quad[9] = { { 0, 0 },
                                        { 0, 1 },
                                        { 0, 2 },
                                        { 0, 3 },
                                        { 1, 0 },
                                        { 1, 1 },
                                        { 1, 2 },
                                        { 1, 3 },
                                        { 2, 0 } };
static const SampleLoc quad_no_edge[9] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 2, 0 } };
static const SampleLoc full_tet[15] = { { 0, 0 },
                                        { 0, 1 },
                                        { 0, 2 },
                                        { 0, 3 },
                                        { 1, 0 },
                                        { 1, 1 },
                                        { 1, 2 },
                                        { 1, 3 },
                                        { 1, 4 },
                                        { 1, 5 },
                                        { 2, 0 },
                                        { 2, 1 },
                                        { 2, 2 },
                                        { 2, 3 },
                                        { 3, 0 } };
static const SampleLoc tet_no_edge[9] = { { 0, 0 },
                                          { 0, 1 },
                                          { 0, 2 },
                                          { 0, 3 },
                                          { 2, 0 },
                                          { 2, 1 },
                                          { 2, 2 },
                                          { 2, 3 },
                                          { 3, 0 } };
static const SampleLoc tet_no_face[11] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 1, 0 },
                                           { 1, 1 },
                                           { 1, 2 },
                                           { 1, 3 },
                                           { 1, 4 },
                                           { 1, 5 },
                                           { 3, 0 } };
static const SampleLoc tet_vtx_vol[5] = { { 0, 0 },
                                          { 0, 1 },
                                          { 0, 2 },
                                          { 0, 3 },
                                          { 3, 0 } };
static const SampleLoc full_pyr[18] = { { 0, 0 },
                                        { 0, 1 },
                                        { 0, 2 },
                                        { 0, 3 },
                                        // Not the apex  {0, 4}
                                        { 1, 0 },
                                        { 1, 1 },
                                        { 1, 2 },
                                        { 1, 3 },
                                        { 1, 4 },
                                        { 1, 5 },
                                        { 1, 6 },
                                        { 1, 7 },
                                        { 2, 0 },
                                        { 2, 1 },
                                        { 2, 2 },
                                        { 2, 3 },
                                        { 2, 4 },
                                        { 3, 0 } };
static const SampleLoc pyr_no_edge[10] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           // Not the apex  {0, 4}
                                           { 2, 0 },
                                           { 2, 1 },
                                           { 2, 2 },
                                           { 2, 3 },
                                           { 2, 4 },
                                           { 3, 0 } };
static const SampleLoc pyr_no_face[13] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           // Not the apex  {0, 4}
                                           { 1, 0 },
                                           { 1, 1 },
                                           { 1, 2 },
                                           { 1, 3 },
                                           { 1, 4 },
                                           { 1, 5 },
                                           { 1, 6 },
                                           { 3, 0 } };
static const SampleLoc pyr_vtx_vol[ 5] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           // Not the apex  {0, 4}
                                           { 3, 0 } };
static const SampleLoc full_wdg[21] = { { 0, 0 },
                                        { 0, 1 },
                                        { 0, 2 },
                                        { 0, 3 },
                                        { 0, 4 },
                                        { 0, 5 },
                                        { 1, 0 },
                                        { 1, 1 },
                                        { 1, 2 },
                                        { 1, 3 },
                                        { 1, 4 },
                                        { 1, 5 },
                                        { 1, 6 },
                                        { 1, 7 },
                                        { 1, 8 },
                                        { 2, 0 },
                                        { 2, 1 },
                                        { 2, 2 },
                                        { 2, 3 },
                                        { 2, 4 },
                                        { 3, 0 } };
static const SampleLoc wdg_no_edge[12] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 0, 4 },
                                           { 0, 5 },
                                           { 2, 0 },
                                           { 2, 1 },
                                           { 2, 2 },
                                           { 2, 3 },
                                           { 2, 4 },
                                           { 3, 0 } };
static const SampleLoc wdg_no_face[16] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 0, 4 },
                                           { 0, 5 },
                                           { 1, 0 },
                                           { 1, 1 },
                                           { 1, 2 },
                                           { 1, 3 },
                                           { 1, 4 },
                                           { 1, 5 },
                                           { 1, 6 },
                                           { 1, 7 },
                                           { 1, 8 },
                                           { 3, 0 } };
static const SampleLoc wdg_vtx_vol[ 7] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 0, 4 },
                                           { 0, 5 },
                                           { 3, 0 } };
static const SampleLoc full_hex[27] = { { 0, 0 },
                                        { 0, 1 },
                                        { 0, 2 },
                                        { 0, 3 },
                                        { 0, 4 },
                                        { 0, 5 },
                                        { 0, 6 },
                                        { 0, 7 },
                                        { 1, 0 },
                                        { 1, 1 },
                                        { 1, 2 },
                                        { 1, 3 },
                                        { 1, 4 },
                                        { 1, 5 },
                                        { 1, 6 },
                                        { 1, 7 },
                                        { 1, 8 },
                                        { 1, 9 },
                                        { 1,10 },
                                        { 1,11 },
                                        { 2, 0 },
                                        { 2, 1 },
                                        { 2, 2 },
                                        { 2, 3 },
                                        { 2, 4 },
                                        { 2, 5 },
                                        { 3, 0 } };
static const SampleLoc hex_no_edge[15] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 0, 4 },
                                           { 0, 5 },
                                           { 0, 6 },
                                           { 0, 7 },
                                           { 2, 0 },
                                           { 2, 1 },
                                           { 2, 2 },
                                           { 2, 3 },
                                           { 2, 4 },
                                           { 2, 5 },
                                           { 3, 0 } };
static const SampleLoc hex_no_face[21] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 0, 4 },
                                           { 0, 5 },
                                           { 0, 6 },
                                           { 0, 7 },
                                           { 1, 0 },
                                           { 1, 1 },
                                           { 1, 2 },
                                           { 1, 3 },
                                           { 1, 4 },
                                           { 1, 5 },
                                           { 1, 6 },
                                           { 1, 7 },
                                           { 1, 8 },
                                           { 1, 9 },
                                           { 1,10 },
                                           { 1,11 },
                                           { 3, 0 } };
static const SampleLoc hex_vtx_vol[ 9] = { { 0, 0 },
                                           { 0, 1 },
                                           { 0, 2 },
                                           { 0, 3 },
                                           { 0, 4 },
                                           { 0, 5 },
                                           { 0, 6 },
                                           { 0, 7 },
                                           { 3, 0 } };
static const SampleLoc invalid[26] = { { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 },
                                       { -1, -1 } };

/* Initialize and return look-up table for location_from_sample_number */
static loc_table_ptr get_location_table() 
{
  static const SampleLoc* location_table[MIXED][16];
  for( unsigned i = 0; i < MIXED; ++i)
    for (unsigned j = 0; j <16; ++j)
      location_table[i][j] = invalid;
  
  location_table[TRIANGLE][1] = full_tri;
  location_table[TRIANGLE][2] = full_tri+3;  // + num_vtx
  location_table[TRIANGLE][3] = full_tri;
  location_table[TRIANGLE][4] = full_tri+6;  // + num_vtx + num_edge
  location_table[TRIANGLE][5] = tri_no_edge;
  location_table[TRIANGLE][6] = full_tri + 3;// + num_vtx
  location_table[TRIANGLE][7] = full_tri;
  
  location_table[QUADRILATERAL][1] = full_quad;
  location_table[QUADRILATERAL][2] = full_quad+4;  // + num_vtx
  location_table[QUADRILATERAL][3] = full_quad;
  location_table[QUADRILATERAL][4] = full_quad+8;  // + num_vtx + num_edge
  location_table[QUADRILATERAL][5] = quad_no_edge;
  location_table[QUADRILATERAL][6] = full_quad + 4;// + num_vtx
  location_table[QUADRILATERAL][7] = full_quad;
  
  location_table[TETRAHEDRON][1] = full_tet;
  location_table[TETRAHEDRON][2] = full_tet+4;   // + num_vtx
  location_table[TETRAHEDRON][3] = full_tet;
  location_table[TETRAHEDRON][4] = full_tet+10;  // + num_vtx + num_edge
  location_table[TETRAHEDRON][5] = tet_no_edge;
  location_table[TETRAHEDRON][6] = full_tet+4;   // + num_vtx
  location_table[TETRAHEDRON][7] = full_tet;
  location_table[TETRAHEDRON][8] = full_tet+14;  // + num_vtx + num_edge + num_face
  location_table[TETRAHEDRON][9] = tet_vtx_vol;
  location_table[TETRAHEDRON][10]= tet_no_face+4;// + num_vtx
  location_table[TETRAHEDRON][11]= tet_no_face;
  location_table[TETRAHEDRON][12]= full_tet+10;  // + num_vtx + num_edge
  location_table[TETRAHEDRON][13]= tet_no_edge;
  location_table[TETRAHEDRON][14]= full_tet+4;   // + num_vtx
  location_table[TETRAHEDRON][15]= full_tet;
  
  location_table[PRISM][1] = full_wdg;
  location_table[PRISM][2] = full_wdg+6;   // + num_vtx
  location_table[PRISM][3] = full_wdg;
  location_table[PRISM][4] = full_wdg+15;  // + num_vtx + num_edge
  location_table[PRISM][5] = wdg_no_edge;
  location_table[PRISM][6] = full_wdg+6;   // + num_vtx
  location_table[PRISM][7] = full_wdg;
  location_table[PRISM][8] = full_wdg+20;  // + num_vtx + num_edge + num_face
  location_table[PRISM][9] = wdg_vtx_vol;
  location_table[PRISM][10]= wdg_no_face+6;// + num_vtx
  location_table[PRISM][11]= wdg_no_face;
  location_table[PRISM][12]= full_wdg+15;  // + num_vtx + num_edge
  location_table[PRISM][13]= wdg_no_edge;
  location_table[PRISM][14]= full_wdg+6;   // + num_vtx
  location_table[PRISM][15]= full_wdg;
   
  location_table[PYRAMID][1] = full_pyr;
  location_table[PYRAMID][2] = full_pyr+4;   // + num_vtx
  location_table[PYRAMID][3] = full_pyr;
  location_table[PYRAMID][4] = full_pyr+12;  // + num_vtx + num_edge
  location_table[PYRAMID][5] = pyr_no_edge;
  location_table[PYRAMID][6] = full_pyr+4;   // + num_vtx
  location_table[PYRAMID][7] = full_pyr;
  location_table[PYRAMID][8] = full_pyr+17;  // + num_vtx + num_edge + num_face
  location_table[PYRAMID][9] = pyr_vtx_vol;
  location_table[PYRAMID][10]= pyr_no_face+4;// + num_vtx
  location_table[PYRAMID][11]= pyr_no_face;
  location_table[PYRAMID][12]= full_pyr+12;  // + num_vtx + num_edge
  location_table[PYRAMID][13]= pyr_no_edge;
  location_table[PYRAMID][14]= full_pyr+4;   // + num_vtx
  location_table[PYRAMID][15]= full_pyr;
   
  location_table[HEXAHEDRON][1] = full_hex;
  location_table[HEXAHEDRON][2] = full_hex+8;   // + num_vtx
  location_table[HEXAHEDRON][3] = full_hex;
  location_table[HEXAHEDRON][4] = full_hex+20;  // + num_vtx + num_edge
  location_table[HEXAHEDRON][5] = hex_no_edge;
  location_table[HEXAHEDRON][6] = full_hex+8;   // + num_vtx
  location_table[HEXAHEDRON][7] = full_hex;
  location_table[HEXAHEDRON][8] = full_hex+26;  // + num_vtx + num_edge + num_face
  location_table[HEXAHEDRON][9] = hex_vtx_vol;
  location_table[HEXAHEDRON][10]= hex_no_face+8;// + num_vtx
  location_table[HEXAHEDRON][11]= hex_no_face;
  location_table[HEXAHEDRON][12]= full_hex+20;  // + num_vtx + num_edge
  location_table[HEXAHEDRON][13]= hex_no_edge;
  location_table[HEXAHEDRON][14]= full_hex+8;   // + num_vtx
  location_table[HEXAHEDRON][15]= full_hex;
 
  return location_table;
}

} // namespace Mesquite
