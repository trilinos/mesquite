/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD: 18-Jun-04 at 12:23:49 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file TargetCalculatorTest.cpp

Unit testing of various TargetCalculator concrete classes. 
*/
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "FeasibleNewton.hpp"
#include "TargetCalculator.hpp"
#include "ConcreteTargetCalculators.hpp"
#include "LPtoPTemplate.hpp"
#include "sRI_DFT.hpp"
#include "MeshImpl.hpp"
#include "InstructionQueue.hpp"

#include "PatchDataInstances.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "MsqFreeVertexIndexIterator.hpp"
#include <list>
#include <iterator>
#include <limits>

using namespace Mesquite;
using std::cout;
using std::endl;
using std::cerr;

class TargetCalculatorTest : public CppUnit::TestFixture,
                             public Mesquite::WTargetCalculator
{
private:
  CPPUNIT_TEST_SUITE(TargetCalculatorTest);
  CPPUNIT_TEST (test_DefaultTargetCalculator);
  CPPUNIT_TEST (test_compute_Lambda);
  CPPUNIT_TEST (test_compute_V);
  CPPUNIT_TEST (test_compute_Q);
  CPPUNIT_TEST (test_compute_Delta);
  CPPUNIT_TEST (test_optimize_vertex_positions_tets);  
  CPPUNIT_TEST (test_optimize_vertex_positions_hexes);  
  CPPUNIT_TEST (test_optimize_vertex_positions_triangles);  
  CPPUNIT_TEST (test_optimize_vertex_positions_quads);  
  CPPUNIT_TEST (test_local_patches);
  CPPUNIT_TEST_SUITE_END();
   
private:

  Matrix3D mG;
  
  PatchData m4Quads;
  PatchData m6Quads;
  PatchData m12Hex;
  PatchData triPatch;
  PatchData tetPatch;

public:
  void setUp()
  {
    MsqPrintError err(cout);

    mG =  "  2.   0.  -1  "
          "  1.   1.  -1  "
          " -1.   0.   2. ";
     
    /* our 2D set up: 4 quads, center vertex outcentered by (0,-0.5)
       7____6____5
       |    |    |
       | 2  |  3 |
       8-_  |  _-4       vertex 1 is at (0,0)
       |  -_0_-  |       vertex 5 is at (2,2)
       | 0  |  1 |
       1----2----3
    */
    create_four_quads_patch(m4Quads, err); CPPUNIT_ASSERT(!err);

    /*! \fn create_six_quads_patch(PatchData &four_quads, MsqError &err)
      our 2D set up: 6 quads, 1 center vertex outcentered by (0,-0.5), the other centered
      7____6____5___11
      |    |    |    |
      | 2  |  3 | 5  |
      8-_  |  _-4---10       vertex 1 is at (0,0)
      |  -_0_-  |    |       vertex 11 is at (3,2)
      | 0  |  1 | 4  |
      1----2----3----9
    */
    create_six_quads_patch_with_domain(m6Quads, err); CPPUNIT_ASSERT(!err);

    /*! \fn create_twelve_hex_patch(PatchData &pd, MsqError &err)
      3D set up: 12 quads, one center vertex outcentered by (0,-0.5),
      the other centered. Vertex 1 is at (0,0,-1). Vertex 35 is at (3,2,1).
     
      7____6____5___11     19___18____17__23     31___30___29___35
      |    |    |    |      |    |    |    |      |    |    |    |
      | 2  |  3 | 5  |      |    |    |    |      | 8  |  9 | 11 |
      8----0----4---10     20-_  |  _16---22     32---24---28---34       
      |    |    |    |      |  -12_-  |    |      |    |    |    |       
      | 0  |  1 | 4  |      |    |    |    |      | 6  |  7 | 10 |
      1----2----3----9     13---14---15---21     25---26---27---33
    */
    create_twelve_hex_patch(m12Hex, err); CPPUNIT_ASSERT(!err);

   /*! \fn create_two_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing two ideal triangles
          / 0 \
         0-----1
          \ 1 /
           \ /
            3
   */
    create_qm_two_tri_patch_with_domain(triPatch,err);CPPUNIT_ASSERT(!err);
    
    create_qm_two_tet_patch(tetPatch,err);CPPUNIT_ASSERT(!err);
    
  }

  void tearDown()
  {
    destroy_patch_with_domain(triPatch);
    destroy_patch_with_domain(m6Quads);
  }
  
public:
  TargetCalculatorTest()
  {}
  
  void test_DefaultTargetCalculator();

  void test_compute_Lambda();

  void test_compute_V();

  void test_compute_Q();
  void test_compute_Delta();

  /*! This test optimizes a mesh with one free vertex and makes sure the vertex moves
      to a given position. The default target matrices and sRI DFT (size and rotation free)
      are used
      \param file_name is the name of the VTK file that contains the mesh.
      \param vtx_index is the index of the vertex that we expect to move.
      \param res the position we expect the vertex to move to.
      \param normal for a surface mesh, this is the normal to the mesh.
      \param point for a surface mesh, this is a point on the surface.
    */
  void test_optimize_vertex_positions(const char* file_name, size_t vtx_index, Vector3D res,
                              Vector3D* normal=0, Vector3D* point=0);

  void test_optimize_vertex_positions_tets();
  
  void test_optimize_vertex_positions_hexes();
  
  void test_optimize_vertex_positions_triangles();
  
  void test_optimize_vertex_positions_quads();
  
  void test_local_patches();
};

  
void TargetCalculatorTest::test_DefaultTargetCalculator()
{
  MsqPrintError err(cout);

  // Creates calculator and compute isotropic target corner matrices.
  ShapeGuides811 iso_calc;
  iso_calc.compute_target_matrices(triPatch, triPatch, err); CPPUNIT_ASSERT(!err);

  TargetMatrix W;

  // checks corner matrices for first triangle, first corner.
  W = triPatch.targetMatrices.get_element_corner_tags( &triPatch, 0, err )[0];
  //W = elems[0].get_tag()->target_matrix(0);
  double fac = pow(2./sqrt(3.), 1./3.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[0][0], fac*1., 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[0][1], fac*.5, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[0][2], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[1][0], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[1][1], fac*sqrt(3.)/2., 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[1][2], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[2][0], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[2][1], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[2][2], fac, 1e-6);

  // checks corner matrices for second triangle, third corner.
  W = triPatch.targetMatrices.get_element_corner_tags( &triPatch, 0, err )[2];
  //W = elems[1].get_tag()->target_matrix(2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[0][0], fac*1., 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[0][1], fac*.5, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[0][2], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[1][0], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[1][1], fac*sqrt(3.)/2., 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[1][2], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[2][0], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[2][1], 0, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(W[2][2], fac, 1e-6);

  // checks there isn't a 4th corner matrix available
  // (normally commented out, since the embedded assert will stop the code). 
  // W = elems[1].get_tag()->target_matrix(3);
}


void TargetCalculatorTest::test_compute_Lambda()
{
  MsqPrintError err(cout);
  
  double Lambda = compute_Lambda(mG, err); CPPUNIT_ASSERT(!err);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4422, Lambda, .0001);

  Lambda = m4Quads.get_average_Lambda_3d(err); CPPUNIT_ASSERT(!err);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.984604, Lambda, .0001); // this test result is not 100% checked
  
  Lambda = m12Hex.get_average_Lambda_3d(err); CPPUNIT_ASSERT(!err);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.994868, Lambda, .0001); // this test result is not 100% checked
  
}

void TargetCalculatorTest::test_compute_V()
{
  MsqPrintError err(cout);
  
  Matrix3D V = compute_V_3D(mG, err); CPPUNIT_ASSERT(!err);

  Matrix3D result = "  0.8165   -0.3651   0.4472 "
                    "  0.4082    0.9129   0.     "
                    " -0.4082    0.1826   0.8944 ";
  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(result[i][j], V[i][j], .0001);
}

void TargetCalculatorTest::test_compute_Q()
{
  MsqPrintError err(cout);
  
  Matrix3D Q = compute_Q_3D(mG, err); CPPUNIT_ASSERT(!err);

  Matrix3D result = "  1.2599    0.5143  -1.0499 "
                    "  0.        1.1502  -0.0939 "
                    "  0.        0.       0.6900 ";
  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(result[i][j], Q[i][j], .0001);
}

void TargetCalculatorTest::test_compute_Delta()
{
  MsqPrintError err(cout);
  
  Matrix3D Delta = compute_Delta_3D(mG, err); CPPUNIT_ASSERT(!err);

  Matrix3D result = "  1.348   0.  0. "
                    "  0.      0.5503  0. "
                    "  0.      0.  1.348 ";

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(result[i][j], Delta[i][j], .0001);
}


/*! This test optimizes a mesh with one free vertex and makes sure the vertex moves
    to a given position. The default target matrices and sRI DFT (size and rotation free)
    are used
    \param file_name is the name of the VTK file that contains the mesh.
    \param vtx_index is the index of the vertex that we expect to move.
    \param res the position we expect the vertex to move to.
    \param normal for a surface mesh, this is the normal to the mesh.
    \param point for a surface mesh, this is a point on the surface.
  */
void TargetCalculatorTest::test_optimize_vertex_positions(
                                             const char* file_name, 
                                             size_t vtx_index, 
                                             Vector3D res,
                                             Vector3D* normal, 
                                             Vector3D* point)
{
  MsqPrintError err(cout);

  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk(file_name, err); CPPUNIT_ASSERT(!err);
  PlanarDomain* domain = 0;
  if (normal) {
    domain = new PlanarDomain(*normal, *point);
  }
  
  sRI_DFT dft;
  LPtoPTemplate obj_func(&dft, 1, err); CPPUNIT_ASSERT(!err);
  FeasibleNewton improver(&obj_func);
  improver.set_patch_type(PatchData::GLOBAL_PATCH, err); CPPUNIT_ASSERT(!err);
  ShapeGuides811 targ_calc;

  InstructionQueue q;
  q.add_target_calculator( &targ_calc, err ); CPPUNIT_ASSERT(!err);
  q.set_master_quality_improver( &improver, err ); CPPUNIT_ASSERT(!err);
  q.run_instructions( mesh, domain, err ); CPPUNIT_ASSERT(!err);

  //improver.set_target_calculator(&targ_calc, err); CPPUNIT_ASSERT(!err);
      
  //improver.loop_over_mesh(mesh_set, err); CPPUNIT_ASSERT(!err);

  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err ); CPPUNIT_ASSERT(!err);
  MsqVertex* vtx = pd.get_vertex_array(err); CPPUNIT_ASSERT(!err);
  int good = vtx[vtx_index].within_tolerance_box(res, 1e-4);
  cout << "vtx[]: " << vtx[vtx_index] << endl;
  CPPUNIT_ASSERT( good==1 );
}

void TargetCalculatorTest::test_optimize_vertex_positions_tets()
{
  Vector3D res(0,0,0);
  size_t vtx_index=0;
  test_optimize_vertex_positions("../../meshFiles/3D/VTK/two_tets_shape.vtk",
                                 vtx_index, res);
}

void TargetCalculatorTest::test_optimize_vertex_positions_hexes()
{
  Vector3D res(0,0,0);
  size_t vtx_index=26;
  test_optimize_vertex_positions("../../meshFiles/3D/VTK/hex_2_with_bound.vtk",
                                 vtx_index, res);
}

void TargetCalculatorTest::test_optimize_vertex_positions_triangles()
{
  Vector3D res(0,0,5);
  size_t vtx_index=8;
  Vector3D normal (0,0,1);
  Vector3D point (0,0,5);
  test_optimize_vertex_positions("../../meshFiles/2D/VTK/square_tri_2.vtk",
                                 vtx_index, res, &normal, &point);
}

void TargetCalculatorTest::test_optimize_vertex_positions_quads()
{
  Vector3D res(0,0,5);
  size_t vtx_index=8;
  Vector3D normal (0,0,1);
  Vector3D point (0,0,5);
  test_optimize_vertex_positions("../../meshFiles/2D/VTK/square_quad_2.vtk",
                                 vtx_index, res, &normal, &point);
}

void TargetCalculatorTest::test_local_patches()
{
  const double EPS = std::numeric_limits<double>::epsilon();
  
    // Could use any file...
  const char* file_name = "../../meshFiles/3D/VTK/large_box_hex_1000.vtk";
  
  MsqPrintError err(cout);
  ShapeGuides811 targ_calc;
  InstructionQueue q;
  q.add_target_calculator( &targ_calc, err ); CPPUNIT_ASSERT(!err);
  
    // Read two copies of the same mesh
  Mesquite::MeshImpl mesh1, mesh2;
  mesh1.read_vtk( file_name, err ); CPPUNIT_ASSERT(!err);
  mesh2.read_vtk( file_name, err ); CPPUNIT_ASSERT(!err);
  
    // Calculate target matrices using global patch for mesh 1
  targ_calc.set_patch_type( PatchData::GLOBAL_PATCH, err );
  CPPUNIT_ASSERT(!err);
  q.run_instructions( &mesh1, err ); 
  CPPUNIT_ASSERT(!err);
  
    // Calulate target matrices using local patch for mesh 2
  targ_calc.set_patch_type( PatchData::ELEMENTS_ON_VERTEX_PATCH, err, 1 );
  CPPUNIT_ASSERT(!err);
  q.run_instructions( &mesh2, err ); 
  CPPUNIT_ASSERT(!err);
  
    // Get global patches for each mesh
  PatchData patch1, patch2;
  patch1.set_mesh( &mesh1 );
  patch2.set_mesh( &mesh2 );
  bool b; size_t junk;
  targ_calc.set_patch_type( PatchData::GLOBAL_PATCH, err ); CPPUNIT_ASSERT(!err);
  b = patch1.get_next_vertex_element_patch( 1, false, junk, err ); CPPUNIT_ASSERT(b && !err);
  b = patch2.get_next_vertex_element_patch( 1, false, junk, err ); CPPUNIT_ASSERT(b && !err);
  CPPUNIT_ASSERT( patch1.num_vertices() == patch2.num_vertices() );
  CPPUNIT_ASSERT( patch1.num_elements() == patch2.num_elements() );
  CPPUNIT_ASSERT( patch1.num_corners() == patch2.num_corners() );
  CPPUNIT_ASSERT( patch1.num_vertices() != 0 );
  CPPUNIT_ASSERT( patch1.num_elements() != 0 );
  CPPUNIT_ASSERT( patch1.num_corners() != 0 );
 
    // Compare target matrices
  const TargetMatrix *m1, *m2;
  unsigned total = 0;
  for (unsigned i = 0; i < patch1.num_elements(); ++i)  // for each element
  {
    m1 = patch1.targetMatrices.get_element_corner_tags( &patch1, i, err );
    CPPUNIT_ASSERT(!err);
    m2 = patch1.targetMatrices.get_element_corner_tags( &patch2, i, err );
    CPPUNIT_ASSERT(!err);
    size_t n = patch1.element_by_index( i ).corner_count();
    total += n;
    for (unsigned j = 0; j < n; ++j)  // for each target matrix in element
    {
        // for each matrix element
      for (unsigned r = 0; r < 3; ++r) 
        for (unsigned c = 0; c < 3; ++c)
        {
          const double d1 = m1[j][r][c];
          const double d2 = m2[j][r][c];
          CPPUNIT_ASSERT_DOUBLES_EQUAL( d1, d2, d1 * EPS );
        }
        // check cK also
      const double cK1 = m1[j].get_cK();
      const double cK2 = m1[j].get_cK();
      CPPUNIT_ASSERT_DOUBLES_EQUAL( cK1, cK2, cK1 * EPS );
    } // for (j)
  } // for(i)
  std::cout << "Compared " << total << " target matrices" << std::endl;
}


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetCalculatorTest, "TargetCalculatorTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetCalculatorTest, "Unit");
