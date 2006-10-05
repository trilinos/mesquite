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
# include <iostream.h>
#else
# include <iostream>
  using std::cout;
  using std::endl;
#endif

#include <string>
#ifdef MSQ_USE_OLD_STD_HEADERS
# include <map.h>
# include <set.h>
# include <vector.h>
#else
# include <map>
# include <set>
# include <vector>
  using std::set;
  using std::map;
  using std::vector;
  using std::string;
#endif


#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "MsqVertex.hpp"
#include "MeshTSTT.hpp"
#include "GeomTSTT.hpp"
#include "cppunit/extensions/HelperMacros.h"


#include "TSTTB.hh"
#include "TSTTM.hh"

using namespace Mesquite;

/* This is the mesh defined in the two arrays below.  This test
 * should work with any mesh composed entirely of triangle
 * elements.  The mesh used for the test can be changed by modifying
 * the two arrays below.
 *                                             
 *            3-------------------2            
 *           / \                 / \           
 *          /   \               /   \          
 *         /     \             /     \         
 *        /       \    (1)    /       \        
 *       /         \         /         \       
 *      /    (2)    \       /    (0)    \      
 *     /             \     /             \     
 *    /               \   /               \    
 *   /                 \ /                 \   
 *  4-------------------0-------------------1  
 *   \                 / \                 /   
 *    \               /   \               /    
 *     \             /     \             /     
 *      \    (3)    /       \    (5)    /      
 *       \         /         \         /       
 *        \       /    (4)    \       /        
 *         \     /             \     /         
 *          \   /               \   /          
 *           \ /                 \ /           
 *            5-------------------6            
 *                                             
 */

extern const double vertexCoords[] = {
   0.00,  0.00,  0.00,  // 0
   2.00,  0.00,  0.00,  // 1
   1.00,  1.73,  0.00,  // 2
  -1.00,  1.73,  0.00,  // 3
  -2.00,  0.00,  0.00,  // 4
  -1.00, -1.73,  0.00,  // 5
   1.00, -1.73,  0.00   // 6
};

extern const int triangleConnectivity[] = {
 0, 1, 2, // 0
 0, 2, 3, // 1
 0, 3, 4, // 2
 0, 4, 5, // 3
 0, 5, 6, // 4
 0, 6, 1  // 5
};



class TSTT_Test : public CppUnit::TestFixture
{
  private:
    
    CPPUNIT_TEST_SUITE( TSTT_Test );
    CPPUNIT_TEST( testVertexIterator );
    CPPUNIT_TEST( testVertexByte );
    CPPUNIT_TEST( testVertexAdjacency );
    CPPUNIT_TEST( testElementConnectivity );
    CPPUNIT_TEST( testElementTopology );
    CPPUNIT_TEST( testIntTag );
    CPPUNIT_TEST( testDoubleTag );
    CPPUNIT_TEST_SUITE_END();

    MeshTSTT* myMesh;
    TSTTM::Mesh myTSTTMesh;
    
    Mesh::VertexHandle vtxIndexToHandle[7];
    Mesh::ElementHandle triIndexToHandle[7];
    map<Mesh::VertexHandle,int> vtxHandleToIndex;
    map<Mesh::ElementHandle,int> triHandleToIndex;

   bool match_triangles( const int* tri1, const Mesh::VertexHandle* tri2_handles );
   bool writeVtkFile( const char* name );

  public:
  
    TSTT_Test()
     : myMesh(0)
      {}
  
    void setUp();
    void tearDown();
    
    void matchVertexCoordinates();
    void matchElementConnectivity();
    void testVertexIterator();
    void testVertexByte();
    void testVertexAdjacency();
    void testElementConnectivity();
    void testElementTopology();
    void testIntTag();
    void testDoubleTag();
    
};

bool TSTT_Test::writeVtkFile( const char* filename )
{
  int i;
  
  FILE* file = fopen( filename, "w" );
  if (!file) {
    perror( filename );
    return false;
  }
  
  fputs( "# vtk DataFile Version 2.0\n"
         "Mesquite Mesh\n"
         "ASCII\n"
         "DATASET UNSTRUCTURED_GRID\n",
         file );
  
  const int num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  fprintf( file, "POINTS %d float\n", num_pts );
  for (i = 0; i < num_pts; ++i)
    fprintf( file, "%f %f %f\n",
                   vertexCoords[3*i    ],
                   vertexCoords[3*i + 1],
                   vertexCoords[3*i + 2] );
  
  int num_tris = sizeof(triangleConnectivity) / (3*sizeof(int));
  fprintf( file, "CELLS %d %d\n", num_tris, 4*num_tris );
  for (i = 0; i < num_tris; ++i)
    fprintf( file, "3 %d %d %d\n",
                   triangleConnectivity[3*i    ],
                   triangleConnectivity[3*i + 1],
                   triangleConnectivity[3*i + 2] );
  
  fprintf( file, "CELL_TYPES %d\n", num_tris );
  for (i = 0; i < num_tris; ++i)
    fprintf( file, "5\n" );
    
  fclose( file );
  return true;
}


void TSTT_Test::setUp()
{
  const char* TMP_FILE_NAME = "hexagon.vtk";
  MsqPrintError err(cout);
  void* root_set = 0;
  
  if (!writeVtkFile(TMP_FILE_NAME)) {
    cout << "Cannot write temporary file: " << TMP_FILE_NAME << endl;
    return;
  }
  
  try {
    TSTTM::Factory factory = TSTTM::Factory::_create();
    myTSTTMesh = factory.newMesh("");
    if (!myTSTTMesh) {
      cout << "Cannot access TSTT mesh implementation." << endl;
      return;
    }
  
    myTSTTMesh.load( 0, TMP_FILE_NAME );
    root_set = myTSTTMesh.getRootSet();
  } 
  catch (...) {
    cout << "TSTT Mesh implementation failed to read file." << endl;
    return;
  }
  
  //remove( TMP_FILE_NAME );
          
  myMesh = MeshTSTT::create( myTSTTMesh, myTSTTMesh.getRootSet(), err );
  if (err || !myMesh)
  {
    delete myMesh;
    myMesh = 0;
    return;
  }
  
  matchVertexCoordinates();
  matchElementConnectivity();
}

void TSTT_Test::tearDown()
{
  delete myMesh;
  myMesh = 0;
}

void TSTT_Test::matchVertexCoordinates()
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );
    // make sure this hasn't been called yet
  CPPUNIT_ASSERT( 0 == vtxHandleToIndex.size() );
    // initialize data
  memset( vtxIndexToHandle, 0, sizeof(vtxIndexToHandle) );
  
    // Get vertex handles
  vector<Mesh::VertexHandle> vertices;
  myMesh->get_all_vertices( vertices, err ); 
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( num_pts, vertices.size() );
  
    // get vertex coordinates
  vector<MsqVertex> coordinates( num_pts );
  myMesh->vertices_get_coordinates( &vertices[0],
                                    &coordinates[0],
                                    num_pts,
                                    err );
  CPPUNIT_ASSERT( !err );
  
    // match vertex coordiantes
  for (size_t i = 0; i < num_pts; ++i)
  {
    Mesquite::Vector3D coord( vertexCoords[3*i], vertexCoords[3*i+1], vertexCoords[3*i+2] );
    size_t j;
    for (j = 0; j < vertices.size(); ++j)
    {
      if (((coordinates[j]) - coord).length() < DBL_EPSILON)
      {
        vtxIndexToHandle[i] = vertices[j];
        CPPUNIT_ASSERT(vtxHandleToIndex.find(vertices[j]) == vtxHandleToIndex.end());
        vtxHandleToIndex[vertices[j]] = i;
        break;
      }
    }
    
    CPPUNIT_ASSERT(j < vertices.size()); // found a match
  }
}

bool TSTT_Test::match_triangles( const int* tri1, 
                                 const Mesh::VertexHandle* tri2_handles )
{
  int tri2[3] = { vtxHandleToIndex[tri2_handles[0]],
                  vtxHandleToIndex[tri2_handles[1]],
                  vtxHandleToIndex[tri2_handles[2]] };
  
  return (tri1[0] == tri2[0] && tri1[1] == tri2[1] && tri1[2] == tri2[2]) ||
         (tri1[0] == tri2[1] && tri1[1] == tri2[2] && tri1[2] == tri2[0]) ||
         (tri1[0] == tri2[2] && tri1[1] == tri2[0] && tri1[2] == tri2[1]);
}

void TSTT_Test::matchElementConnectivity()
{
  MsqPrintError err(cout);
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );
    // make sure this hasn't been called yet
  CPPUNIT_ASSERT( 0 == triHandleToIndex.size() );
    // initialize data
  memset( triIndexToHandle, 0, sizeof(triIndexToHandle) );
  
  vector<Mesh::VertexHandle> vertices;
  vector<Mesh::ElementHandle> elements;
  vector<size_t> offsets;
  myMesh->get_all_elements( elements, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( elements.size() == num_tri );
  myMesh->elements_get_attached_vertices( &elements[0],
                                          elements.size(),
                                          vertices,
                                          offsets,
                                          err );
  CPPUNIT_ASSERT(!err);
                                          
    // Make sure all are triangles
  size_t i;
  for (i = 0; i < elements.size(); ++i)
    CPPUNIT_ASSERT( offsets[i] + 3 == offsets[i+1] );
  
    // Match triangles
  for (size_t i = 0; i < num_tri; ++i)
  {
    size_t j;
    for (j = 0; j < elements.size(); ++j)
    {
      Mesh::VertexHandle verts[3] = {
        vertices[offsets[j]  ],
        vertices[offsets[j]+1],
        vertices[offsets[j]+2] };
      
      if (match_triangles( triangleConnectivity + 3*i, verts ))
      {
        triIndexToHandle[i] = elements[j];
        CPPUNIT_ASSERT(triHandleToIndex.find(elements[j]) == triHandleToIndex.end());
        triHandleToIndex[elements[j]] = i;
        break;
      }
    }
    
    CPPUNIT_ASSERT(j < elements.size()); // found a match
  }
}


void TSTT_Test::testVertexIterator() 
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // mark each vertex as it is encountered
  msq_std::vector<int> marks(num_pts);
  memset( &marks[0], 0, num_pts * sizeof(int) );
  
    // iterate over vertices
  unsigned count = 0;
  VertexIterator* iter = myMesh->vertex_iterator(err);
  CPPUNIT_ASSERT(!err && iter);
  while (!iter->is_at_end())
  {
    Mesh::VertexHandle handle = iter->operator*();
    iter->operator++();
    
    map<Mesh::VertexHandle,int>::iterator f = vtxHandleToIndex.find(handle);
    CPPUNIT_ASSERT( f != vtxHandleToIndex.end() );
    
    unsigned index = f->second;
    CPPUNIT_ASSERT(index < num_pts);
    CPPUNIT_ASSERT(marks[index] == 0);
    marks[index] = 1;
    ++count;
  }
  CPPUNIT_ASSERT(count == num_pts);
}  

void TSTT_Test::testVertexByte()
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // set each one individually
  unsigned char bytes[num_pts];
  size_t i;
  for (i = 0; i < num_pts; ++i)
  {
    bytes[i] = (unsigned char)(rand() % 128);
    myMesh->vertex_set_byte( vtxIndexToHandle[i], bytes[i], err );
    CPPUNIT_ASSERT( !err );
  }
  for (i = 0; i < num_pts; ++i)
  {
    unsigned char byte = 0;
    myMesh->vertex_get_byte( vtxIndexToHandle[i], &byte, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT( byte == bytes[i] );
  }
  
    // set all at once
  for (i = 0; i < num_pts; ++i)
    bytes[i] = (unsigned char)(rand() % 128);
  myMesh->vertices_set_byte( vtxIndexToHandle, bytes, num_pts, err );
  CPPUNIT_ASSERT( !err );
  unsigned char bytes2[num_pts];
  myMesh->vertices_get_byte( vtxIndexToHandle, bytes2, num_pts, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( !memcmp( bytes, bytes2, num_pts ) );
}  
  

void TSTT_Test::testVertexAdjacency()
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // construct adjacency list from input data to compare against
  msq_std::set<int> adjset[num_pts];
  size_t i;
  for (i = 0; i < num_tri; ++i)
    for (size_t j = 3*i; j < 3*i+3; ++j)
      adjset[triangleConnectivity[j]].insert(i);
  
  msq_std::vector<Mesh::ElementHandle> elements;
  msq_std::vector<size_t> offsets;
  myMesh->vertices_get_attached_elements( vtxIndexToHandle, num_pts,
                                          elements, offsets, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(offsets.size() == num_pts + 1);
  
    // compare connectivity for each vertex
  for (i = 0; i < num_pts; ++i)
  {
    size_t count = offsets[i+1] - offsets[i];
    CPPUNIT_ASSERT(adjset[i].size() == count);
    Mesh::ElementHandle* elems = &elements[offsets[i]];
    
    for (size_t j = 0; j < count; ++j)
    {
        // Get element index from handle
      Mesh::ElementHandle handle = elems[j];
      msq_std::map<Mesh::ElementHandle,int>::iterator idx = triHandleToIndex.find(handle);
      CPPUNIT_ASSERT(idx != triHandleToIndex.end());

        // Find element index in vertex adjacency set
      msq_std::set<int>::iterator iter = adjset[i].find(idx->second);
      CPPUNIT_ASSERT(iter != adjset[i].end());
      adjset[i].erase(iter);
    }
      // Make sure we got all the adjacent elements
    CPPUNIT_ASSERT(adjset[i].size() == 0);
  }
}


void TSTT_Test::testElementConnectivity()
{
  MsqPrintError err(cout);
  //const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // get element connectivity list
  msq_std::vector<Mesh::VertexHandle> vertices;
  msq_std::vector<size_t> offsets;
  myMesh->elements_get_attached_vertices( triIndexToHandle, num_tri,
                                          vertices, offsets, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(offsets.size() == num_tri + 1);
  CPPUNIT_ASSERT(vertices.size() == 3*num_tri);
  
    // check each element's connectivity
  Mesh::VertexHandle elem_vertices[3];
  for (size_t i = 0; i < num_tri; ++i)
  {
      // check that connectivity list contains three vertices
    CPPUNIT_ASSERT(offsets[i] + 3 == offsets[i+1]);
    
      // get list of vertex indices from connectivity data
    for (size_t j = 0; j < 3; j++)
    {
      size_t offset = offsets[i] + j;
      CPPUNIT_ASSERT( offset < vertices.size() );
      elem_vertices[j] = vertices[offset];
    }
    
      // compare connectivity
    CPPUNIT_ASSERT( match_triangles( triangleConnectivity + 3*i, elem_vertices ) );
  }
}
    
    
void TSTT_Test::testElementTopology()
{
  MsqPrintError err(cout);
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

  EntityTopology topo[num_tri];
  myMesh->elements_get_topologies( triIndexToHandle, topo, num_tri, err );
  CPPUNIT_ASSERT(!err);
  for (size_t i = 0; i < num_tri; ++i)
    CPPUNIT_ASSERT( topo[i] == Mesquite::TRIANGLE );
}


void TSTT_Test::testIntTag()
{
  const char* tagname = "TEST_TEST_INT_TAG";
  MsqPrintError err(cout);
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // create a tag
  TagHandle tag = myMesh->tag_create( tagname, Mesh::INT, 2, NULL, err );
  CPPUNIT_ASSERT( !err );
  
    // get the tag
  TagHandle tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( tag == tag2 );
  
    // check tag metadata
  string name;
  Mesh::TagType type;
  unsigned length;
  myMesh->tag_properties( tag, name, type, length, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( name == tagname );
  CPPUNIT_ASSERT( type == Mesh::INT );
  CPPUNIT_ASSERT( length == 2 );
  
    // set the tag on all triangles
  std::vector<int> data1(2*num_tri), data2(2*num_tri);
  std::vector<int>::iterator iter1, iter2;
  for (iter1 = data1.begin(); iter1 != data1.end(); ++iter1)
    *iter1 = rand();
  myMesh->tag_set_element_data( tag, num_tri, triIndexToHandle, &data1[0], err );
  CPPUNIT_ASSERT( !err );
  
    // get tag data from all triangles and compare
  myMesh->tag_get_element_data( tag, num_tri, triIndexToHandle, &data2[0], err );
  CPPUNIT_ASSERT( !err );
  for (iter1 = data1.begin(), iter2 = data2.begin(); 
       iter1 != data1.end(); ++iter1, ++iter2)
    CPPUNIT_ASSERT( *iter1 == *iter2 );
    
    // destroy the tag
  myMesh->tag_delete( tag, err );
  CPPUNIT_ASSERT(!err);
  tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT(err.error_code() == MsqError::TAG_NOT_FOUND);
  err.clear();
}

  
void TSTT_Test::testDoubleTag()
{
  const char* tagname = "TEST_TEST_DOUBLE_TAG";
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // create a tag
  TagHandle tag = myMesh->tag_create( tagname, Mesh::DOUBLE, 1, NULL, err );
  CPPUNIT_ASSERT( !err );
  
    // get the tag
  TagHandle tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( tag == tag2 );
  
    // check tag metadata
  string name;
  Mesh::TagType type;
  unsigned length;
  myMesh->tag_properties( tag, name, type, length, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( name == tagname );
  CPPUNIT_ASSERT( type == Mesh::DOUBLE );
  CPPUNIT_ASSERT( length == 1 );
  
    // set the tag on all vertices
  std::vector<double> data1(num_pts), data2(num_pts);
  std::vector<double>::iterator iter1, iter2;
  for (iter1 = data1.begin(); iter1 != data1.end(); ++iter1)
    *iter1 = sqrt(abs(rand()));
  myMesh->tag_set_vertex_data( tag, num_pts, vtxIndexToHandle, &data1[0], err );
  CPPUNIT_ASSERT( !err );
  
    // get tag data from all vertices and compare
  myMesh->tag_get_vertex_data( tag, num_pts, vtxIndexToHandle, &data2[0], err );
  CPPUNIT_ASSERT( !err );
  for (iter1 = data1.begin(), iter2 = data2.begin(); 
       iter1 != data1.end(); ++iter1, ++iter2)
    CPPUNIT_ASSERT( *iter1 == *iter2 );
    
    // destroy the tag
  myMesh->tag_delete( tag, err );
  CPPUNIT_ASSERT(!err);
  tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT(err.error_code() == MsqError::TAG_NOT_FOUND);
  err.clear();
}

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TSTT_Test, "TSTT_Test");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TSTT_Test, "Unit");
