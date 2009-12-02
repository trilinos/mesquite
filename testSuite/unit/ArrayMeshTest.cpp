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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ArrayMeshTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ArrayMesh.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "UnitUtil.hpp"

using namespace Mesquite;

class ArrayMeshTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE( ArrayMeshTest );
  
  CPPUNIT_TEST( test_get_geometric_dimension );
  CPPUNIT_TEST( test_get_all_elements );
  CPPUNIT_TEST( test_get_all_vertices );
  CPPUNIT_TEST( test_get_all_vertices_one_based );
  CPPUNIT_TEST( test_vertex_iterator );
  CPPUNIT_TEST( test_vertex_iterator_one_based );
  CPPUNIT_TEST( test_element_iterator );
  CPPUNIT_TEST( test_vertices_get_fixed_flag );
  CPPUNIT_TEST( test_vertices_get_coordinates );
  CPPUNIT_TEST( test_vertices_get_coordinates_two_d );
  CPPUNIT_TEST( test_vertices_get_coordinates_one_based );
  CPPUNIT_TEST( test_vertex_set_coordinates );
  CPPUNIT_TEST( test_vertex_set_coordinates_two_d );
  CPPUNIT_TEST( test_vertex_set_coordinates_one_based );
  CPPUNIT_TEST( test_vertex_set_byte );
  CPPUNIT_TEST( test_vertices_set_byte );
  CPPUNIT_TEST( test_vertices_get_attached_elements );
  CPPUNIT_TEST( test_vertices_get_attached_elements_one_based );
  CPPUNIT_TEST( test_elements_get_attached_vertices );
  CPPUNIT_TEST( test_elements_get_topologies );

  CPPUNIT_TEST( test_vertices_get_attached_elements_mixed );
  CPPUNIT_TEST( test_vertices_get_attached_elements_mixed_one_based );
  CPPUNIT_TEST( test_elements_get_attached_vertices_mixed );
  CPPUNIT_TEST( test_elements_get_attached_vertices_mixed_one_based );
  CPPUNIT_TEST( test_elements_get_topologies_mixed );
  
  CPPUNIT_TEST_SUITE_END();

  ArrayMesh *zeroBased3D, *oneBased3D, *zeroBased2D, *oneBased2D;
  ArrayMesh *mixedZeroBased, *mixedOneBased;
  double *zeroBased3Dcoords, *oneBased3Dcoords, *zeroBased2Dcoords, *oneBased2Dcoords;
  double *mixedZeroBasedCoords, *mixedOneBasedCoords;
  
public:

  void setUp();
  void tearDown();

  void test_get_geometric_dimension ();
  void test_get_all_elements ();
  void test_get_all_vertices ();
  void test_get_all_vertices_one_based ();
  void test_vertex_iterator ();
  void test_vertex_iterator_one_based ();
  void test_element_iterator ();
  void test_vertices_get_fixed_flag ();
  void test_vertices_get_coordinates ();
  void test_vertices_get_coordinates_two_d ();
  void test_vertices_get_coordinates_one_based ();
  void test_vertex_set_coordinates ();
  void test_vertex_set_coordinates_two_d ();
  void test_vertex_set_coordinates_one_based ();
  void test_vertex_set_byte ();
  void test_vertices_set_byte ();
  void test_vertices_get_attached_elements ();
  void test_vertices_get_attached_elements_one_based ();
  void test_elements_get_attached_vertices ();
  void test_elements_get_topologies ();

  void test_vertices_get_attached_elements_mixed ();
  void test_vertices_get_attached_elements_mixed_one_based ();
  void test_elements_get_attached_vertices_mixed ();
  void test_elements_get_attached_vertices_mixed_one_based ();
  void test_elements_get_topologies_mixed ();
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ArrayMeshTest, "ArrayMeshTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ArrayMeshTest, "Unit");

/* Mesh:

   0---------1----------2
   |         |          |
   |    0    |    1     |
   |         |          |
   |         |          |
   3---------4----------5
   |         |          |
   Y    2    |    3     |
   ^         |          |
   |         |          |
   6-->X-----7----------8

 Z = -2
   
*/
const double coords_2d[] = { 0, 2,
                             1, 2,
                             2, 2,
                             0, 1,
                             1, 1,
                             2, 1,
                             0, 0,
                             1, 0,
                             2, 0 };
const double coords_3d[] = { 0, 2, -2,
                             1, 2, -2,
                             2, 2, -2,
                             0, 1, -2,
                             1, 1, -2,
                             2, 1, -2,
                             0, 0, -2,
                             1, 0, -2,
                             2, 0, -2};
const unsigned long conn_zero_based[] = { 3, 4, 1, 0,
                                          4, 5, 2, 1,
                                          6, 7, 4, 3, 
                                          7, 8, 5, 4 };
const unsigned long conn_one_based[] = {  4, 5, 2, 1,
                                          5, 6, 3, 2,
                                          7, 8, 5, 4, 
                                          8, 9, 6, 5 };

const int fixed[] = { 1, 1, 1,
                      1, 0, 1,
                      1, 1, 1 };

/* Mixed mesh
 *
 *   0-------1-------2-------3
 *   |       | \  2  |       |
 *   |   0   |   \   |   3   |
 *   |       | 1   \ |       |
 *   7-------6-------5-------4
 */
const double mixed_coords[] = { 0.0, 1.0, 0.0,
                                1.0, 1.0, 0.0,
                                2.0, 1.0, 0.0,
                                3.0, 1.0, 0.0,
                                3.0, 0.0, 0.0,
                                2.0, 0.0, 0.0,
                                1.0, 0.0, 0.0,
                                0.0, 0.0, 0.0 };
const unsigned long mixed_conn_zero[] = { 7, 6, 1, 0,
                                          6, 5, 1,
                                          2, 1, 5,
                                          5, 4, 3, 2 };
const unsigned long mixed_conn_one[]  = { 8, 7, 2, 1,
                                          7, 6, 2,
                                          3, 2, 6,
                                          6, 5, 4, 3 };
const EntityTopology mixed_types[] = { QUADRILATERAL,
                                       TRIANGLE,
                                       TRIANGLE,
                                       QUADRILATERAL };
const unsigned long conn_offsets_zero[] = { 0, 4, 7, 10, 14 };


void ArrayMeshTest::setUp()
{
     zeroBased3Dcoords = new double[27];
      oneBased3Dcoords = new double[27];
     zeroBased2Dcoords = new double[18];
      oneBased2Dcoords = new double[18];
  mixedZeroBasedCoords = new double[24];
   mixedOneBasedCoords = new double[24];
  memcpy(    zeroBased3Dcoords, coords_3d,    27*sizeof(double) );
  memcpy(     oneBased3Dcoords, coords_3d,    27*sizeof(double) );
  memcpy(    zeroBased2Dcoords, coords_2d,    18*sizeof(double) );
  memcpy(     oneBased2Dcoords, coords_2d,    18*sizeof(double) );
  memcpy( mixedZeroBasedCoords, mixed_coords, 24*sizeof(double) );
  memcpy(  mixedOneBasedCoords, mixed_coords, 24*sizeof(double) );
     zeroBased3D = new ArrayMesh( 3, 9,    zeroBased3Dcoords, fixed, 4, QUADRILATERAL, conn_zero_based,                    false );
      oneBased3D = new ArrayMesh( 3, 9,     oneBased3Dcoords, fixed, 4, QUADRILATERAL, conn_one_based,                     true  );
     zeroBased2D = new ArrayMesh( 2, 9,    zeroBased2Dcoords, fixed, 4, QUADRILATERAL, conn_zero_based,                    false );
      oneBased2D = new ArrayMesh( 2, 9,     oneBased2Dcoords, fixed, 4, QUADRILATERAL, conn_one_based,                     true  );
  mixedZeroBased = new ArrayMesh( 3, 8, mixedZeroBasedCoords, fixed, 4, mixed_types,   mixed_conn_zero, conn_offsets_zero, false );
   mixedOneBased = new ArrayMesh( 3, 8,  mixedOneBasedCoords, fixed, 4, mixed_types,   mixed_conn_one,  NULL,              true  );
}

void ArrayMeshTest::tearDown()
{
  delete    zeroBased3D;
  delete     oneBased3D;
  delete    zeroBased2D;
  delete     oneBased2D;
  delete mixedZeroBased;
  delete  mixedOneBased;
  delete []    zeroBased3Dcoords;
  delete []     oneBased3Dcoords;
  delete []    zeroBased2Dcoords;
  delete []     oneBased2Dcoords;
  delete [] mixedZeroBasedCoords;
  delete []  mixedOneBasedCoords;
}

void ArrayMeshTest::test_get_geometric_dimension()
{
  MsqPrintError err(std::cerr);
  CPPUNIT_ASSERT_EQUAL( 2, zeroBased2D->get_geometric_dimension( err ) );
  CPPUNIT_ASSERT_EQUAL( 2,  oneBased2D->get_geometric_dimension( err ) );
  CPPUNIT_ASSERT_EQUAL( 3, zeroBased3D->get_geometric_dimension( err ) );
  CPPUNIT_ASSERT_EQUAL( 3,  oneBased3D->get_geometric_dimension( err ) );
}

void ArrayMeshTest::test_get_all_elements()
{
  MsqPrintError err(std::cerr);
  std::vector<Mesh::ElementHandle> list;
  zeroBased3D->get_all_elements( list, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( (size_t)4, list.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)list[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)list[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)list[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)list[3] );
}

void ArrayMeshTest::test_get_all_vertices()
{
  MsqPrintError err(std::cerr);
  std::vector<Mesh::VertexHandle> list;
  zeroBased3D->get_all_vertices( list, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( (size_t)9, list.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)list[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)list[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)list[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)list[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t)4, (size_t)list[4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)5, (size_t)list[5] );
  CPPUNIT_ASSERT_EQUAL( (size_t)6, (size_t)list[6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)7, (size_t)list[7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)8, (size_t)list[8] );
}

void ArrayMeshTest::test_get_all_vertices_one_based()
{
  MsqPrintError err(std::cerr);
  std::vector<Mesh::VertexHandle> list;
  oneBased3D->get_all_vertices( list, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( (size_t)9, list.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)list[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)list[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)list[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)4, (size_t)list[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t)5, (size_t)list[4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)6, (size_t)list[5] );
  CPPUNIT_ASSERT_EQUAL( (size_t)7, (size_t)list[6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)8, (size_t)list[7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)9, (size_t)list[8] );
}

void ArrayMeshTest::test_vertex_iterator()
{
  MsqPrintError err(std::cerr);
  VertexIterator* iter = zeroBased3D->vertex_iterator( err );
  CPPUNIT_ASSERT( !err );
  std::auto_ptr<VertexIterator> deleter(iter);
  
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)4, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)5, (size_t)iter->operator*() );
   iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)6, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)7, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)8, (size_t)iter->operator*() );
  
  iter->restart();
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)iter->operator*() );
}

void ArrayMeshTest::test_vertex_iterator_one_based()
{
  MsqPrintError err(std::cerr);
  VertexIterator* iter = oneBased3D->vertex_iterator( err );
  CPPUNIT_ASSERT( !err );
  std::auto_ptr<VertexIterator> deleter(iter);
  
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)4, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)5, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)6, (size_t)iter->operator*() );
   iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)7, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)8, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)9, (size_t)iter->operator*() );
  
  iter->restart();
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)iter->operator*() );
}

void ArrayMeshTest::test_element_iterator()
{
  MsqPrintError err(std::cerr);
  ElementIterator* iter = zeroBased3D->element_iterator( err );
  CPPUNIT_ASSERT( !err );
  std::auto_ptr<ElementIterator> deleter(iter);
  
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(!iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)iter->operator*() );
  iter->operator++();
  CPPUNIT_ASSERT(iter->is_at_end());
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)iter->operator*() );
  
  iter->restart();
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)iter->operator*() );
}

void ArrayMeshTest::test_vertices_get_fixed_flag()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 4 };
  bool flags[10];
  zeroBased3D->vertices_get_fixed_flag( (const Mesh::VertexHandle*)verts,
                                        flags, 10, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(  true, flags[0] );
  CPPUNIT_ASSERT_EQUAL(  true, flags[1] );
  CPPUNIT_ASSERT_EQUAL(  true, flags[2] );
  CPPUNIT_ASSERT_EQUAL(  true, flags[3] );
  CPPUNIT_ASSERT_EQUAL( false, flags[4] );
  CPPUNIT_ASSERT_EQUAL(  true, flags[5] );
  CPPUNIT_ASSERT_EQUAL(  true, flags[6] );
  CPPUNIT_ASSERT_EQUAL(  true, flags[7] );
  CPPUNIT_ASSERT_EQUAL(  true, flags[8] );
  CPPUNIT_ASSERT_EQUAL( false, flags[9] );
}

void ArrayMeshTest::test_vertices_get_coordinates()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 4 };
  MsqVertex coords[10];
  zeroBased3D->vertices_get_coordinates( (const Mesh::VertexHandle*)verts,
                                          coords, 10, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 0), coords[0], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 3), coords[1], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 6), coords[2], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 9), coords[3], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+12), coords[4], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+15), coords[5], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+18), coords[6], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+21), coords[7], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+24), coords[8], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+12), coords[4], DBL_EPSILON );
}

void ArrayMeshTest::test_vertices_get_coordinates_two_d()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 4 };
  MsqVertex coords[10];
  zeroBased2D->vertices_get_coordinates( (const Mesh::VertexHandle*)verts,
                                          coords, 10, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[ 0],coords_2d[ 1],0), coords[0], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[ 2],coords_2d[ 3],0), coords[1], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[ 4],coords_2d[ 5],0), coords[2], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[ 6],coords_2d[ 7],0), coords[3], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[ 8],coords_2d[ 9],0), coords[4], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[10],coords_2d[11],0), coords[5], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[12],coords_2d[13],0), coords[6], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[14],coords_2d[15],0), coords[7], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[16],coords_2d[17],0), coords[8], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_2d[ 8],coords_2d[ 9],0), coords[9], DBL_EPSILON );
}

void ArrayMeshTest::test_vertices_get_coordinates_one_based()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 5 };
  MsqVertex coords[10];
  oneBased3D->vertices_get_coordinates( (const Mesh::VertexHandle*)verts,
                                          coords, 10, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 0), coords[0], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 3), coords[1], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 6), coords[2], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+ 9), coords[3], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+12), coords[4], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+15), coords[5], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+18), coords[6], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+21), coords[7], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+24), coords[8], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(coords_3d+12), coords[4], DBL_EPSILON );
}

void ArrayMeshTest::test_vertex_set_coordinates()
{
  MsqPrintError err(std::cerr);
  Mesh::VertexHandle vert = (Mesh::VertexHandle)4;
  Vector3D new_pos( 5, 5, 5 );
  MsqVertex pos;
  zeroBased3D->vertex_set_coordinates( vert, new_pos, err );
  CPPUNIT_ASSERT(!err);
  zeroBased3D->vertices_get_coordinates( &vert, &pos, 1, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( new_pos, pos, DBL_EPSILON );
}

void ArrayMeshTest::test_vertex_set_coordinates_two_d()
{
  MsqPrintError err(std::cerr);
  Mesh::VertexHandle vert = (Mesh::VertexHandle)4;
  Vector3D new_pos( 5, 5, 5 );
  MsqVertex pos;
  zeroBased2D->vertex_set_coordinates( vert, new_pos, err );
  CPPUNIT_ASSERT(!err);
  zeroBased2D->vertices_get_coordinates( &vert, &pos, 1, err );
  CPPUNIT_ASSERT(!err);
  new_pos[2] = 0;
  CPPUNIT_ASSERT_VECTORS_EQUAL( new_pos, pos, DBL_EPSILON );
}

void ArrayMeshTest::test_vertex_set_coordinates_one_based()
{
  MsqPrintError err(std::cerr);
  Mesh::VertexHandle vert = (Mesh::VertexHandle)5;
  Vector3D new_pos( 5, 5, 5 );
  MsqVertex pos;
  oneBased3D->vertex_set_coordinates( vert, new_pos, err );
  CPPUNIT_ASSERT(!err);
  oneBased3D->vertices_get_coordinates( &vert, &pos, 1, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( new_pos, pos, DBL_EPSILON );
}

void ArrayMeshTest::test_vertex_set_byte()
{
  MsqPrintError err(std::cerr);
  unsigned char byte = 'A';
  Mesh::VertexHandle vert = (Mesh::VertexHandle)5;
  zeroBased3D->vertex_set_byte( vert, byte, err );
  CPPUNIT_ASSERT(!err);
  byte = '0';
  zeroBased3D->vertex_get_byte( vert, &byte, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (unsigned char)'A', byte );
}

void ArrayMeshTest::test_vertices_set_byte()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 4 };
  const unsigned char bytes[] = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'E' };
  unsigned char b[10];
  memset( b, 0, 10 );
  zeroBased3D->vertices_set_byte( (const Mesh::VertexHandle*)verts, bytes, 10, err );
  CPPUNIT_ASSERT(!err);
  zeroBased3D->vertices_get_byte( (const Mesh::VertexHandle*)verts, b, 10, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( bytes[0], b[0] );
  CPPUNIT_ASSERT_EQUAL( bytes[1], b[1] );
  CPPUNIT_ASSERT_EQUAL( bytes[2], b[2] );
  CPPUNIT_ASSERT_EQUAL( bytes[3], b[3] );
  CPPUNIT_ASSERT_EQUAL( bytes[4], b[4] );
  CPPUNIT_ASSERT_EQUAL( bytes[5], b[5] );
  CPPUNIT_ASSERT_EQUAL( bytes[6], b[6] );
  CPPUNIT_ASSERT_EQUAL( bytes[7], b[7] );
  CPPUNIT_ASSERT_EQUAL( bytes[8], b[8] );
}


void ArrayMeshTest::test_vertices_get_attached_elements()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
  std::vector<Mesh::ElementHandle> elems;
  std::vector<size_t> offsets;
  zeroBased3D->vertices_get_attached_elements( (const Mesh::VertexHandle*)verts,
                                               9, elems, offsets, err );
  CPPUNIT_ASSERT(!err);

  CPPUNIT_ASSERT_EQUAL( (size_t)16, elems.size() );
  CPPUNIT_ASSERT( offsets.size() == 9 || offsets.size() == 10 );
  if (offsets.size() == 10)
    CPPUNIT_ASSERT_EQUAL( size_t(16), offsets[9] );

    // vert 0
  CPPUNIT_ASSERT_EQUAL( (size_t)0, offsets[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 0] );
    // vert 1
  CPPUNIT_ASSERT_EQUAL( (size_t)1, offsets[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 2] );
    // vert 2
  CPPUNIT_ASSERT_EQUAL( (size_t)3, offsets[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 3] );
    // vert 3
  CPPUNIT_ASSERT_EQUAL( (size_t)4, offsets[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 5] );
    // vert 4
  CPPUNIT_ASSERT_EQUAL( (size_t)6, offsets[4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 8] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 9] );
    // vert 5
  CPPUNIT_ASSERT_EQUAL( (size_t)10, offsets[5] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[10] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[11] );
    // vert 6
  CPPUNIT_ASSERT_EQUAL( (size_t)12, offsets[6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[12] );
    // vert 7
  CPPUNIT_ASSERT_EQUAL( (size_t)13, offsets[7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[13] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[14] );
    // vert 8
  CPPUNIT_ASSERT_EQUAL( (size_t)15, offsets[8] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[15] );
}

void ArrayMeshTest::test_vertices_get_attached_elements_one_based()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  std::vector<Mesh::ElementHandle> elems;
  std::vector<size_t> offsets;
  oneBased3D->vertices_get_attached_elements( (const Mesh::VertexHandle*)verts,
                                               9, elems, offsets, err );
  CPPUNIT_ASSERT(!err);

  CPPUNIT_ASSERT_EQUAL( (size_t)16, elems.size() );
  CPPUNIT_ASSERT( offsets.size() == 9 || offsets.size() == 10 );
  if (offsets.size() == 10)
    CPPUNIT_ASSERT_EQUAL( size_t(16), offsets[9] );

    // vert 0
  CPPUNIT_ASSERT_EQUAL( (size_t)0, offsets[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 0] );
    // vert 1
  CPPUNIT_ASSERT_EQUAL( (size_t)1, offsets[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 2] );
    // vert 2
  CPPUNIT_ASSERT_EQUAL( (size_t)3, offsets[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 3] );
    // vert 3
  CPPUNIT_ASSERT_EQUAL( (size_t)4, offsets[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 5] );
    // vert 4
  CPPUNIT_ASSERT_EQUAL( (size_t)6, offsets[4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 8] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 9] );
    // vert 5
  CPPUNIT_ASSERT_EQUAL( (size_t)10, offsets[5] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[10] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[11] );
    // vert 6
  CPPUNIT_ASSERT_EQUAL( (size_t)12, offsets[6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[12] );
    // vert 7
  CPPUNIT_ASSERT_EQUAL( (size_t)13, offsets[7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[13] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[14] );
    // vert 8
  CPPUNIT_ASSERT_EQUAL( (size_t)15, offsets[8] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[15] );
}

void ArrayMeshTest::test_elements_get_attached_vertices()
{
  MsqPrintError err(std::cerr);
  const size_t elems[] = { 3, 2, 1, 0 };
  std::vector<Mesh::VertexHandle> verts;
  std::vector<size_t> offsets;
  zeroBased3D->elements_get_attached_vertices( (const Mesh::ElementHandle*)elems,
                                               4, verts, offsets, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (size_t)16,  verts.size() );
  CPPUNIT_ASSERT( offsets.size() == 4 || offsets.size() == 5 );
  if (offsets.size() == 5)
    CPPUNIT_ASSERT_EQUAL( (size_t)16, offsets[4] );
  
    // elem 3
  CPPUNIT_ASSERT_EQUAL( (size_t) 0, offsets[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 7, (size_t)verts[ 0] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 8, (size_t)verts[ 1] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 5, (size_t)verts[ 2] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, (size_t)verts[ 3] );
  
    // elem 2
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, offsets[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 6, (size_t)verts[ 4] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 7, (size_t)verts[ 5] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, (size_t)verts[ 6] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 3, (size_t)verts[ 7] );
   
    // elem 1
  CPPUNIT_ASSERT_EQUAL( (size_t) 8, offsets[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, (size_t)verts[ 8] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 5, (size_t)verts[ 9] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 2, (size_t)verts[10] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 1, (size_t)verts[11] );
   
    // elem 0
  CPPUNIT_ASSERT_EQUAL( (size_t)12, offsets[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 3, (size_t)verts[12] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, (size_t)verts[13] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 1, (size_t)verts[14] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 0, (size_t)verts[15] );
} 


void ArrayMeshTest::test_elements_get_topologies()
{
  MsqPrintError err(std::cerr);
  const size_t elems[] = { 3, 2, 1, 0, 1, 2, 3, 0 };
  const size_t num_elem = sizeof(elems)/sizeof(elems[0]);
  EntityTopology topo[num_elem];
  memset( topo, 0, sizeof(topo) );
  zeroBased3D->elements_get_topologies( (const Mesh::ElementHandle*)elems, 
                                         topo, num_elem, err );
  CPPUNIT_ASSERT(!err);
  for (size_t i = 0; i < num_elem; ++i)
    CPPUNIT_ASSERT_EQUAL( (int)QUADRILATERAL, (int)topo[i] );
}


void ArrayMeshTest::test_vertices_get_attached_elements_mixed()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  std::vector<Mesh::ElementHandle> elems;
  std::vector<size_t> offsets;
  mixedZeroBased->vertices_get_attached_elements( (const Mesh::VertexHandle*)verts,
                                               8, elems, offsets, err );
  CPPUNIT_ASSERT(!err);

  CPPUNIT_ASSERT_EQUAL( (size_t)14, elems.size() );
  CPPUNIT_ASSERT( offsets.size() == 8 || offsets.size() == 9 );
  if (offsets.size() == 9)
    CPPUNIT_ASSERT_EQUAL( elems.size(), offsets[8] );

    // vert 0
  CPPUNIT_ASSERT_EQUAL( (size_t)0, offsets[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 0] );
    // vert 1
  CPPUNIT_ASSERT_EQUAL( (size_t)1, offsets[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 3] );
    // vert 2
  CPPUNIT_ASSERT_EQUAL( (size_t)4, offsets[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 5] );
    // vert 3
  CPPUNIT_ASSERT_EQUAL( (size_t)6, offsets[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 6] );
    // vert 4
  CPPUNIT_ASSERT_EQUAL( (size_t)7, offsets[4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 7] );
    // vert 5
  CPPUNIT_ASSERT_EQUAL( (size_t)8, offsets[5] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[8] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[9] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[10] );
    // vert 6
  CPPUNIT_ASSERT_EQUAL( (size_t)11, offsets[6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[11] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[12] );
    // vert 7
  CPPUNIT_ASSERT_EQUAL( (size_t)13, offsets[7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[13] );
}

void ArrayMeshTest::test_vertices_get_attached_elements_mixed_one_based()
{
  MsqPrintError err(std::cerr);
  const size_t verts[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  std::vector<Mesh::ElementHandle> elems;
  std::vector<size_t> offsets;
  mixedOneBased->vertices_get_attached_elements( (const Mesh::VertexHandle*)verts,
                                               8, elems, offsets, err );
  CPPUNIT_ASSERT(!err);

  CPPUNIT_ASSERT_EQUAL( (size_t)14, elems.size() );
  CPPUNIT_ASSERT( offsets.size() == 8 || offsets.size() == 9 );
  if (offsets.size() == 9)
    CPPUNIT_ASSERT_EQUAL( elems.size(), offsets[8] );

    // vert 0
  CPPUNIT_ASSERT_EQUAL( (size_t)0, offsets[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 0] );
    // vert 1
  CPPUNIT_ASSERT_EQUAL( (size_t)1, offsets[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[ 1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[ 2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 3] );
    // vert 2
  CPPUNIT_ASSERT_EQUAL( (size_t)4, offsets[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[ 4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 5] );
    // vert 3
  CPPUNIT_ASSERT_EQUAL( (size_t)6, offsets[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 6] );
    // vert 4
  CPPUNIT_ASSERT_EQUAL( (size_t)7, offsets[4] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[ 7] );
    // vert 5
  CPPUNIT_ASSERT_EQUAL( (size_t)8, offsets[5] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[8] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, (size_t)elems[9] );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, (size_t)elems[10] );
    // vert 6
  CPPUNIT_ASSERT_EQUAL( (size_t)11, offsets[6] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[11] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, (size_t)elems[12] );
    // vert 7
  CPPUNIT_ASSERT_EQUAL( (size_t)13, offsets[7] );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, (size_t)elems[13] );
}

void ArrayMeshTest::test_elements_get_attached_vertices_mixed()
{
  MsqPrintError err(std::cerr);
  const size_t elems[] = { 3, 2, 1, 0 };
  std::vector<Mesh::VertexHandle> verts;
  std::vector<size_t> offsets;
  mixedZeroBased->elements_get_attached_vertices( (const Mesh::ElementHandle*)elems,
                                               4, verts, offsets, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (size_t)14,  verts.size() );
  CPPUNIT_ASSERT( offsets.size() == 4 || offsets.size() == 5 );
  if (offsets.size() == 5)
    CPPUNIT_ASSERT_EQUAL( verts.size(), offsets[4] );
  
    // elem 3
  CPPUNIT_ASSERT_EQUAL( (size_t) 0, offsets[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 5, (size_t)verts[ 0] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, (size_t)verts[ 1] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 3, (size_t)verts[ 2] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 2, (size_t)verts[ 3] );
  
    // elem 2
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, offsets[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 2, (size_t)verts[ 4] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 1, (size_t)verts[ 5] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 5, (size_t)verts[ 6] );
   
    // elem 1
  CPPUNIT_ASSERT_EQUAL( (size_t) 7, offsets[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 6, (size_t)verts[ 7] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 5, (size_t)verts[ 8] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 1, (size_t)verts[ 9] );
   
    // elem 0
  CPPUNIT_ASSERT_EQUAL( (size_t)10, offsets[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 7, (size_t)verts[10] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 6, (size_t)verts[11] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 1, (size_t)verts[12] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 0, (size_t)verts[13] );
} 

void ArrayMeshTest::test_elements_get_attached_vertices_mixed_one_based()
{
  MsqPrintError err(std::cerr);
  const size_t elems[] = { 3, 2, 1, 0 };
  std::vector<Mesh::VertexHandle> verts;
  std::vector<size_t> offsets;
  mixedOneBased->elements_get_attached_vertices( (const Mesh::ElementHandle*)elems,
                                               4, verts, offsets, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (size_t)14,  verts.size() );
  CPPUNIT_ASSERT( offsets.size() == 4 || offsets.size() == 5 );
  if (offsets.size() == 5)
    CPPUNIT_ASSERT_EQUAL( verts.size(), offsets[4] );
  
    // elem 3
  CPPUNIT_ASSERT_EQUAL( (size_t) 0, offsets[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 6, (size_t)verts[ 0] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 5, (size_t)verts[ 1] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, (size_t)verts[ 2] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 3, (size_t)verts[ 3] );
  
    // elem 2
  CPPUNIT_ASSERT_EQUAL( (size_t) 4, offsets[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 3, (size_t)verts[ 4] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 2, (size_t)verts[ 5] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 6, (size_t)verts[ 6] );
   
    // elem 1
  CPPUNIT_ASSERT_EQUAL( (size_t) 7, offsets[2] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 7, (size_t)verts[ 7] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 6, (size_t)verts[ 8] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 2, (size_t)verts[ 9] );
   
    // elem 0
  CPPUNIT_ASSERT_EQUAL( (size_t)10, offsets[3] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 8, (size_t)verts[10] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 7, (size_t)verts[11] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 2, (size_t)verts[12] );
  CPPUNIT_ASSERT_EQUAL( (size_t) 1, (size_t)verts[13] );
} 


void ArrayMeshTest::test_elements_get_topologies_mixed()
{
  MsqPrintError err(std::cerr);
  const size_t elems[] = { 3, 2, 1, 0, 1, 2, 3 };
  const size_t num_elem = sizeof(elems)/sizeof(elems[0]);
  EntityTopology topo[num_elem];
  memset( topo, 0, sizeof(topo) );
  mixedZeroBased->elements_get_topologies( (const Mesh::ElementHandle*)elems, 
                                         topo, num_elem, err );
  CPPUNIT_ASSERT(!err);
  for (size_t i = 0; i < num_elem; ++i)
    CPPUNIT_ASSERT_EQUAL( (int)mixed_types[elems[i]], (int)topo[i] );
}
