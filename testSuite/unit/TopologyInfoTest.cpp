/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else
#include <iostream>
using std::cout;
#endif
#include "cppunit/extensions/HelperMacros.h"

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "TopologyInfo.hpp"

using namespace Mesquite;

class TopologyInfoTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(TopologyInfoTest);

   CPPUNIT_TEST (tri);
   CPPUNIT_TEST (tri3);
   CPPUNIT_TEST (tri4);
   CPPUNIT_TEST (tri6);
   CPPUNIT_TEST (tri7);
   
   CPPUNIT_TEST (quad);
   CPPUNIT_TEST (quad4);
   CPPUNIT_TEST (quad5);
   CPPUNIT_TEST (quad8);
   CPPUNIT_TEST (quad9);
   
   CPPUNIT_TEST (tet);
   CPPUNIT_TEST (tet4);
   CPPUNIT_TEST (tet5);
   CPPUNIT_TEST (tet8);
   CPPUNIT_TEST (tet9);
   CPPUNIT_TEST (tet10);
   CPPUNIT_TEST (tet11);
   CPPUNIT_TEST (tet14);
   CPPUNIT_TEST (tet15);

   CPPUNIT_TEST (hex);
   CPPUNIT_TEST (hex8);
   CPPUNIT_TEST (hex9);
   CPPUNIT_TEST (hex14);
   CPPUNIT_TEST (hex15);
   CPPUNIT_TEST (hex20);
   CPPUNIT_TEST (hex21);
   CPPUNIT_TEST (hex26);
   CPPUNIT_TEST (hex27);
 
   CPPUNIT_TEST (pyramid);
   CPPUNIT_TEST (pyramid5);
   CPPUNIT_TEST (pyramid13);
   
   CPPUNIT_TEST (wedge);
   CPPUNIT_TEST (wedge6);
   CPPUNIT_TEST (wedge15);

   CPPUNIT_TEST (polygon);
   CPPUNIT_TEST (polyhedron);
   
   CPPUNIT_TEST (bad_type);

   CPPUNIT_TEST_SUITE_END();

public:

  void setUp() {}
  
  void tearDown() {}
  
  TopologyInfoTest() {}
  
  bool compare_edge( const unsigned* a, const unsigned* b );
  
  bool compare_face( unsigned len, const unsigned* a, const unsigned* b );
  
  bool compare_vol( unsigned len, const unsigned* a, const unsigned* b );
  
  void test_face_elem( EntityTopology topo, 
                       unsigned num_nodes,
                       unsigned num_sides );

  
  void test_vol_elem( EntityTopology topo, 
                      unsigned num_nodes,
                      unsigned num_verts,
                      unsigned num_edges,
                      unsigned num_faces );
  
  void test_poly( EntityTopology topo );
  
  void tri();
  void tri3();
  void tri4();
  void tri6();
  void tri7();
  
  void quad();
  void quad4();
  void quad5();
  void quad8();
  void quad9();
  
  void tet();
  void tet4();
  void tet5();
  void tet8();
  void tet9();
  void tet10();
  void tet11();
  void tet14();
  void tet15();

  void hex();
  void hex8();
  void hex9();
  void hex14();
  void hex15();
  void hex20();
  void hex21();
  void hex26();
  void hex27();
  
  void pyramid();
  void pyramid5();
  void pyramid13();
  
  void wedge();
  void wedge6();
  void wedge15();
  
  void polygon();
  void polyhedron();

  void bad_type();
  
};


  
bool TopologyInfoTest::compare_edge( const unsigned* a, const unsigned* b )
  { return (a[0] == b[0] && a[1] == b[1]) ||
           (a[0] == b[1] && a[1] == b[0]); }

bool TopologyInfoTest::compare_face( unsigned len, const unsigned* a, const unsigned* b )
{
  unsigned i, j;
  for (i = 0; i < len; ++i)
  {
    for (j = 0; j < len; ++j)
    {
      if (a[j] != b[(i+j)%len])
        break;
    }
    if (j == len)
      return true;
  }
  return false;
}

bool TopologyInfoTest::compare_vol( unsigned len, const unsigned* a, const unsigned* b )
{
  for (unsigned i = 0; i < len; ++i)
    if (a[i] != b[i])
      return false;
  return true;
}

void TopologyInfoTest::test_face_elem( EntityTopology topo, 
                     unsigned num_nodes,
                     unsigned num_sides )
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  unsigned index = 0;
  TopologyInfo::higher_order( topo, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!vol);

  unsigned nodes = num_sides + edge * num_sides + face;
  CPPUNIT_ASSERT( num_nodes == nodes );

  unsigned side, dim;
  for (index = 0; index < num_sides; ++index)
  {
    TopologyInfo::side_number( topo, num_nodes, index, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 0 && side == index);
  }

  if (edge)
  {
    for (unsigned s = 0; s < num_sides; ++s)
    {
      TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(dim == 1 && side == s);
    }
  }
  if (face)
  {
    TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 2 && side == 0);
  }
}


void TopologyInfoTest::test_vol_elem( EntityTopology topo, 
                    unsigned num_nodes,
                    unsigned num_verts,
                    unsigned num_edges,
                    unsigned num_faces )
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  unsigned index = 0;
  TopologyInfo::higher_order( topo, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);

  unsigned nodes = num_verts + edge * num_edges + face * num_faces + vol;
  CPPUNIT_ASSERT( num_nodes == nodes );

  unsigned side, dim;
  for (index = 0; index < num_verts; ++index)
  {
    TopologyInfo::side_number( topo, num_nodes, index, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 0 && side == index);
  }

  if (edge)
  {
    for (unsigned s = 0; s < num_edges; ++s)
    {
      TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(dim == 1 && side == s);
    }
  }
  if (face)
  {
    for (unsigned s = 0; s < num_faces; ++s)
    {
      TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(dim == 2 && side == s);
    }
  }
  if (vol)
  {
    TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 3 && side == 0);
  }
}





void TopologyInfoTest::tri()
{
  MsqPrintError err(cout);

  CPPUNIT_ASSERT (2 == TopologyInfo::dimension( TRIANGLE ));
  CPPUNIT_ASSERT (3 == TopologyInfo::adjacent( TRIANGLE, 1 ));
  CPPUNIT_ASSERT (3 == TopologyInfo::adjacent( TRIANGLE, 0 ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( TRIANGLE, 2 ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( TRIANGLE, 3 ));

  CPPUNIT_ASSERT (3 == TopologyInfo::sides( TRIANGLE ));
  CPPUNIT_ASSERT (3 == TopologyInfo::corners( TRIANGLE ));
  CPPUNIT_ASSERT (3 == TopologyInfo::edges( TRIANGLE ));
  CPPUNIT_ASSERT (0 == TopologyInfo::faces( TRIANGLE ));

  const unsigned num_edges = 3;
  const unsigned* side;
  const unsigned edges[num_edges][2] = { {0, 1}, {1, 2}, {2, 0} };
  unsigned count;
  const unsigned face[] = { 0, 1, 2 };

  for (unsigned i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( TRIANGLE, i, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( TRIANGLE, 1, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  }

  side = TopologyInfo::side_vertices( TRIANGLE, 2, 0, count, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(3 == count);
  CPPUNIT_ASSERT( compare_face( 3, side, face ) );
}

void TopologyInfoTest::tri3()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 3;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}

void TopologyInfoTest::tri4()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 4;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}

void TopologyInfoTest::tri6()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 6;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}

void TopologyInfoTest::tri7()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 7;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}


void TopologyInfoTest::quad()
{
  MsqPrintError err(cout);

  CPPUNIT_ASSERT (2 == TopologyInfo::dimension( QUADRILATERAL ));
  CPPUNIT_ASSERT (4 == TopologyInfo::adjacent( QUADRILATERAL, 1 ));
  CPPUNIT_ASSERT (4 == TopologyInfo::adjacent( QUADRILATERAL, 0 ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( QUADRILATERAL, 2 ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( QUADRILATERAL, 3 ));

  CPPUNIT_ASSERT (4 == TopologyInfo::sides( QUADRILATERAL ));
  CPPUNIT_ASSERT (4 == TopologyInfo::corners( QUADRILATERAL ));
  CPPUNIT_ASSERT (4 == TopologyInfo::edges( QUADRILATERAL ));
  CPPUNIT_ASSERT (0 == TopologyInfo::faces( QUADRILATERAL ));

  const unsigned num_edges = 4;
  const unsigned* side;
  const unsigned edges[num_edges][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
  unsigned count;
  const unsigned face[] = { 0, 1, 2, 3 };

  for (unsigned i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( QUADRILATERAL, i, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( QUADRILATERAL, 1, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  }

  side = TopologyInfo::side_vertices( QUADRILATERAL, 2, 0, count, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( 4 == count );
  CPPUNIT_ASSERT( compare_face( 4, side, face ) );
}

void TopologyInfoTest::quad4()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 4;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}

void TopologyInfoTest::quad5()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 5;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}

void TopologyInfoTest::quad8()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 8;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}

void TopologyInfoTest::quad9()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 9;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}


void TopologyInfoTest::tet()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 4;
  const unsigned num_edges = 6;
  const unsigned num_faces = 4;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( TETRAHEDRON ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( TETRAHEDRON, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( TETRAHEDRON, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( TETRAHEDRON, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( TETRAHEDRON, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( TETRAHEDRON ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( TETRAHEDRON ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( TETRAHEDRON ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( TETRAHEDRON ));

  const unsigned* side;
  unsigned i, count;
  const unsigned vert_per_face = 3;
  unsigned edges[num_edges][2] = { {0, 1}, {1, 2}, {2, 0},
                                   {0, 3}, {1, 3}, {2, 3} };
  unsigned faces[num_faces][vert_per_face] = { { 0, 1, 3 }, { 1, 2, 3 }, 
                                               { 2, 0, 3 }, { 2, 1, 0 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( TETRAHEDRON, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( TETRAHEDRON, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  for (i = 0; i < num_faces; ++i)
  {
    side = TopologyInfo::face_vertices( TETRAHEDRON, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );

    side = TopologyInfo::side_vertices( TETRAHEDRON, 2, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );
  } 

  side = TopologyInfo::side_vertices( TETRAHEDRON, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest::tet4()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 4;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet5()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 5;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet8()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 8;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet9()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 9;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet10()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 10;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet11()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 11;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet14()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 14;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet15()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 15;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}


void TopologyInfoTest::hex()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 8;
  const unsigned num_edges = 12;
  const unsigned num_faces = 6;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( HEXAHEDRON ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( HEXAHEDRON, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( HEXAHEDRON, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( HEXAHEDRON, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( HEXAHEDRON, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( HEXAHEDRON ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( HEXAHEDRON ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( HEXAHEDRON ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( HEXAHEDRON ));

  const unsigned* side;
  unsigned i, count;
  const unsigned vert_per_face = 4;
  unsigned edges[num_edges][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                                   { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
                                   { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } };
  unsigned faces[num_faces][vert_per_face] = { { 0, 1, 5, 4 }, 
                                               { 1, 2, 6, 5 },
                                               { 2, 3, 7, 6 },
                                               { 3, 0, 4, 7 },
                                               { 3, 2, 1, 0 },
                                               { 4, 5, 6, 7 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( HEXAHEDRON, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( HEXAHEDRON, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  for (i = 0; i < num_faces; ++i)
  {
    side = TopologyInfo::face_vertices( HEXAHEDRON, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );

    side = TopologyInfo::side_vertices( HEXAHEDRON, 2, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );
  } 

  side = TopologyInfo::side_vertices( HEXAHEDRON, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest:: hex8()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 8;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest:: hex9()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 9;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest:: hex14()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 14;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex15()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 15;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest:: hex20()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 20;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex21()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 21;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex26()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 26;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex27()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 27;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest::pyramid()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 5;
  const unsigned num_edges = 8;
  const unsigned num_faces = 5;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( PYRAMID ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( PYRAMID, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( PYRAMID, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( PYRAMID, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( PYRAMID, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( PYRAMID ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( PYRAMID ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( PYRAMID ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( PYRAMID ));

  const unsigned* side;
  unsigned i, count;
  const unsigned num_tri_faces = 4;
  const unsigned num_quad_faces = 1;
  const bool tri_before_quad = true;
  CPPUNIT_ASSERT( num_tri_faces + num_quad_faces == num_faces );
  unsigned edges[num_edges][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                                   { 0, 4 }, { 1, 4 }, { 2, 4 }, { 3, 4 } };
  unsigned tris[num_tri_faces][3] = { { 0, 1, 4 }, { 1, 2, 4 }, 
                                       { 2, 3, 4 }, { 3, 0, 4 } };
  unsigned quads[num_quad_faces][4] = { { 3, 2, 1, 0 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( PYRAMID, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( PYRAMID, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  const unsigned tri_off = num_quad_faces * !tri_before_quad;
  for (i = 0; i < num_tri_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PYRAMID, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );

    side = TopologyInfo::side_vertices( PYRAMID, 2, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );
  } 

  const unsigned quad_off = num_tri_faces * tri_before_quad;
  for (i = 0; i < num_quad_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PYRAMID, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );

    side = TopologyInfo::side_vertices( PYRAMID, 2, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );
  } 

  side = TopologyInfo::side_vertices( PYRAMID, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest::pyramid5()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 5;

  TopologyInfo::higher_order( PYRAMID, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PYRAMID, num_nodes, 5, 8, 5 );
}

void TopologyInfoTest::pyramid13()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 13;

  TopologyInfo::higher_order( PYRAMID, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PYRAMID, num_nodes, 5, 8, 5 );
}

void TopologyInfoTest::wedge()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 6;
  const unsigned num_edges = 9;
  const unsigned num_faces = 5;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( PRISM ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( PRISM, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( PRISM, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( PRISM, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( PRISM, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( PRISM ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( PRISM ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( PRISM ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( PRISM ));

  const unsigned* side;
  unsigned i, count;
  const unsigned num_tri_faces = 2;
  const unsigned num_quad_faces = 3;
  const bool tri_before_quad = false;
  CPPUNIT_ASSERT( num_tri_faces + num_quad_faces == num_faces );
  unsigned edges[num_edges][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, 
                                   { 0, 3 }, { 1, 4 }, { 2, 5 },
                                   { 3, 4 }, { 4, 5 }, { 5, 3 } };
  unsigned tris[num_tri_faces][3] = { { 2, 1, 0 }, { 3, 4, 5 } };
  unsigned quads[num_quad_faces][4] = { { 0, 1, 4, 3 },
                                        { 1, 2, 5, 4 },
                                        { 2, 0, 3, 5 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( PRISM, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( PRISM, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  const unsigned tri_off = num_quad_faces * !tri_before_quad;
  for (i = 0; i < num_tri_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PRISM, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );

    side = TopologyInfo::side_vertices( PRISM, 2, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );
  } 

  const unsigned quad_off = num_tri_faces * tri_before_quad;
  for (i = 0; i < num_quad_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PRISM, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );

    side = TopologyInfo::side_vertices( PRISM, 2, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );
  } 

  side = TopologyInfo::side_vertices( PRISM, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest::wedge6()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 6;

  TopologyInfo::higher_order( PRISM, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PRISM, num_nodes, 6, 9, 5 );
}

void TopologyInfoTest::wedge15()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 15;

  TopologyInfo::higher_order( PRISM, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PRISM, num_nodes, 6, 9, 5 );
}

void TopologyInfoTest::test_poly( EntityTopology topo )
{
  CPPUNIT_ASSERT( TopologyInfo::adjacent(topo, 1) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::sides(topo) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::corners(topo) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::edges(topo) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::faces(topo) == 0 );
}


void TopologyInfoTest::polygon()
{
  CPPUNIT_ASSERT( TopologyInfo::dimension( POLYGON ) == 2 );
  test_poly( POLYGON );
}

void TopologyInfoTest::polyhedron()
{
  CPPUNIT_ASSERT( TopologyInfo::dimension( POLYHEDRON ) == 3 );
  test_poly( POLYHEDRON );
}


void TopologyInfoTest::bad_type()
{
  Mesquite::MsqError err;
  EntityTopology bad_types[] = { (EntityTopology)0,
                                 (EntityTopology)1,
                                 MIXED,
                                 (EntityTopology)(MIXED + 1)
                               };

  for (unsigned i = 0; i < (sizeof(bad_types)/sizeof(EntityTopology)); ++i)
  {
    CPPUNIT_ASSERT( TopologyInfo::dimension(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::adjacent(bad_types[i], 1) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::sides(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::corners(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::edges(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::faces(bad_types[i]) == 0 );

    bool a,b, c;
    TopologyInfo::higher_order( bad_types[i], 20, a, b, c, err );
    CPPUNIT_ASSERT(err);

    const unsigned* ptr;
    ptr = TopologyInfo::edge_vertices( bad_types[i], 0, err );
    CPPUNIT_ASSERT(err);

    unsigned count;
    ptr = TopologyInfo::face_vertices( bad_types[i], 0, count, err );
    CPPUNIT_ASSERT(err);

    for (unsigned j = 0; j < 4; ++j)
    {
      ptr = TopologyInfo::side_vertices( bad_types[i], j, 0, count, err );
      CPPUNIT_ASSERT(err);
    }

    unsigned dim, side;
    for (unsigned idx = 0; idx < 20; ++idx)
    {
      TopologyInfo::side_number( bad_types[i], idx, idx/2, dim, side, err );
      CPPUNIT_ASSERT(err);
    }
  }
}


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TopologyInfoTest, "TopologyInfoTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TopologyInfoTest, "Unit");


