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
  
  void test_DefaultTargetCalculator()
  {
    MsqPrintError err(cout);

    // Creates calculator and compute isotropic target corner matrices.
    ShapeGuides811 iso_calc;
    iso_calc.compute_target_matrices(triPatch, err); CPPUNIT_ASSERT(!err);

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


  void test_compute_Lambda()
  {
    MsqPrintError err(cout);
    
    double Lambda = compute_Lambda(mG, err); CPPUNIT_ASSERT(!err);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4422, Lambda, .0001);

    Lambda = m4Quads.get_average_Lambda_3d(err); CPPUNIT_ASSERT(!err);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.984604, Lambda, .0001); // this test result is not 100% checked
    
    Lambda = m12Hex.get_average_Lambda_3d(err); CPPUNIT_ASSERT(!err);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.994868, Lambda, .0001); // this test result is not 100% checked
    
  }

  void test_compute_V()
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

  void test_compute_Q()
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

  void test_compute_Delta()
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
  void test_optimize_vertex_positions(const char* file_name, size_t vtx_index, Vector3D res,
                              Vector3D* normal=0, Vector3D* point=0)
  {
    MsqPrintError err(cout);

    Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
    mesh->read_vtk(file_name, err); CPPUNIT_ASSERT(!err);
    MeshSet mesh_set;
    mesh_set.add_mesh(mesh, err); CPPUNIT_ASSERT(!err);
    if (normal) {
      PlanarDomain* domain = new PlanarDomain(*normal, *point);
      mesh_set.set_domain_constraint(domain, err); CPPUNIT_ASSERT(!err);
    }
    
    sRI_DFT dft;
    LPtoPTemplate obj_func(&dft, 1, err); CPPUNIT_ASSERT(!err);
    FeasibleNewton improver(&obj_func);
    improver.set_patch_type(PatchData::GLOBAL_PATCH, err); CPPUNIT_ASSERT(!err);
    ShapeGuides811 targ_calc;

    InstructionQueue q;
    q.add_target_calculator( &targ_calc, err ); CPPUNIT_ASSERT(!err);
    q.set_master_quality_improver( &improver, err ); CPPUNIT_ASSERT(!err);
    q.run_instructions( mesh_set, err ); CPPUNIT_ASSERT(!err);

    //improver.set_target_calculator(&targ_calc, err); CPPUNIT_ASSERT(!err);
        
    //improver.loop_over_mesh(mesh_set, err); CPPUNIT_ASSERT(!err);

    PatchData pd;
    mesh_set.get_next_patch(pd, &improver, err); CPPUNIT_ASSERT(!err);
    MsqVertex* vtx = pd.get_vertex_array(err); CPPUNIT_ASSERT(!err);
    int good = vtx[vtx_index].within_tolerance_box(res, 1e-4);
    cout << "vtx[]: " << vtx[vtx_index] << endl;
    CPPUNIT_ASSERT( good==1 );
  }

  void test_optimize_vertex_positions_tets()
  {
    Vector3D res(0,0,0);
    size_t vtx_index=0;
    test_optimize_vertex_positions("../../meshFiles/3D/VTK/two_tets_shape.vtk",
                                   vtx_index, res);
  }
  
  void test_optimize_vertex_positions_hexes()
  {
    Vector3D res(0,0,0);
    size_t vtx_index=26;
    test_optimize_vertex_positions("../../meshFiles/3D/VTK/hex_2_with_bound.vtk",
                                   vtx_index, res);
  }
  
  void test_optimize_vertex_positions_triangles()
  {
    Vector3D res(0,0,5);
    size_t vtx_index=8;
    Vector3D normal (0,0,1);
    Vector3D point (0,0,5);
    test_optimize_vertex_positions("../../meshFiles/2D/VTK/square_tri_2.vtk",
                                   vtx_index, res, &normal, &point);
  }
  
  void test_optimize_vertex_positions_quads()
  {
    Vector3D res(0,0,5);
    size_t vtx_index=8;
    Vector3D normal (0,0,1);
    Vector3D point (0,0,5);
    test_optimize_vertex_positions("../../meshFiles/2D/VTK/square_quad_2.vtk",
                                   vtx_index, res, &normal, &point);
  }
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetCalculatorTest, "TargetCalculatorTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetCalculatorTest, "Unit");
