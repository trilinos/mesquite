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
/*! \file DistanceFromTargetTest.cpp

Unit testing of various DistanceFromTarget (base *and* concrete classes) functions. 

\author Thomas Leurent
\date 2004-09-31
 */
// DESCRIP-END.
//

#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "PatchDataInstances.hpp"

#include "ConcreteTargetCalculators.hpp"
#include "DistanceFromTarget.hpp"
#include "I_DFT.hpp"
#include "sI_DFT.hpp"
#include "RI_DFT.hpp"
#include "sRI_DFT.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

#include <math.h>

using namespace Mesquite;

using std::cout;
using std::cerr;
using std::endl;

class DistanceFromTargetTest : public CppUnit::TestFixture, Mesquite::DistanceFromTarget
{
private:
  CPPUNIT_TEST_SUITE(DistanceFromTargetTest);
  CPPUNIT_TEST (test_tri_with_default_target_matrix);
  CPPUNIT_TEST (test_quad_with_default_target_matrix);
  CPPUNIT_TEST (test_tet_with_default_target_matrix);
  CPPUNIT_TEST (test_compute_T_matrices);
  CPPUNIT_TEST_SUITE_END();
  
private:

  PatchData oneTri;
  PatchData triPatch;
  PatchData quadPatch;
  PatchData tetPatch;
  PatchData hexPatch;
  PatchData invertedTri;
  PatchData invertedTet;
  PatchData idealTri;
  PatchData idealTet;
    //Tol used for double comparisons
  double qualTol;
  int pF;//PRINT_FLAG
public:
  void setUp()
  {
      //pF=1;//PRINT_FLAG IS ON
    pF=0;//PRINT_FLAG IS OFF
    MsqError err;
    
    qualTol = 1.e-12;
    
     /* Our triangular patch is made of two tris.  tri_1 is a perfect
        equilateral (the ideal for most metrics).  tri_2 is an arbitrary
        triangle.
     */
    create_qm_two_tri_patch_with_domain(triPatch, err);MSQ_CHKERR(err);
    
     /* Our quad patch is made of two quads.  quad_1 is a perfect
        square (the ideal for most metrics).  quad_2 is an arbitrary
        quad.
     */
    create_qm_two_quad_patch_with_domain(quadPatch,err);MSQ_CHKERR(err);
    
    create_one_tet_patch(tetPatch,err);MSQ_CHKERR(err);
    
     /* Our hex patch is made of two hexes.  hex_1 is a perfect
        unit cube (the ideal for most metrics).  hex_2 is an arbitrary
        hex.
     */
     create_qm_two_hex_patch(hexPatch,err);MSQ_CHKERR(err);

       //'ideal' inverted tet
     create_one_inverted_tet_patch(invertedTet, err);MSQ_CHKERR(err);
       //ideal tri
     create_one_tri_patch(idealTri, err);MSQ_CHKERR(err);
       //ideal tet
     create_one_tet_patch(idealTet, err);MSQ_CHKERR(err);

     // Creates a patch of one triangle, with specific coordinates, without domain.
     oneTri.set_num_vertices(3);
     oneTri.set_num_elements(1);
     double coords[3];
     MsqVertex* vert_array = oneTri.get_vertex_array(err);
     coords[0] = 1; coords[1] = -1; coords[2] = -1;
     vert_array[0] = coords;
     coords[0] = 3; coords[1] = 0; coords[2] = -2;
     vert_array[1] = coords;
     coords[0] = 2; coords[1] = -3; coords[2] = 0;
     vert_array[2] = coords;
     size_t indices[3] = { 0, 1, 2 };
     MsqMeshEntity& tri = oneTri.element_by_index(0);
     tri.set_element_type(Mesquite::TRIANGLE);
     memcpy(tri.get_modifiable_vertex_index_array(),
            indices,
            3*sizeof(size_t));

  }

  void tearDown()
  {
    destroy_patch_with_domain(triPatch);
    destroy_patch_with_domain(quadPatch);
  }
  
public:
  DistanceFromTargetTest()
    {}

#undef __FUNC__
#define __FUNC__ "DistanceFromTargetTest::test_tri_with_default_target_matrix"
   void test_tri_with_default_target_matrix()
   {
     MsqError err;
     Matrix3D tri_m3d, quad_m3d, tet_m3d, hex_m3d; 
        
     TargetCalculator::initialize_default_target_matrices(tri_m3d, quad_m3d, tet_m3d, hex_m3d);

     idealTri.allocate_corner_matrices(err); MSQ_CHKERR(err);
     MsqMeshEntity* elem = idealTri.get_element_array(err);
     MsqTag* tag = elem[0].get_tag();
     tag->target_matrix(0) = tri_m3d;
     tag->target_matrix(1) = tri_m3d;
     tag->target_matrix(2) = tri_m3d;

     double coords[3];
     MsqVertex* vert_array = idealTri.get_vertex_array(err);
     double fac = pow(2./sqrt(3.), 1./3.);
     coords[0] = 0; coords[1] = 0; coords[2] = 0;
     vert_array[0] = coords;
     coords[0] = fac; coords[1] = 0; coords[2] = 0;
     vert_array[1] = coords;
     coords[0] = fac*1./2.; coords[1] = fac*sqrt(3.)/2.; coords[2] = 0;
     vert_array[2] = coords;

     Matrix3D T[MSQ_MAX_NUM_VERT_PER_ENT];
     double c_k[MSQ_MAX_NUM_VERT_PER_ENT];
     this->compute_T_matrices(elem[0], idealTri, T, 3, c_k, err); MSQ_CHKERR(err);

//    cout << T[0] << endl;

     Matrix3D Id = " 1 0 0 "
                   " 0 1 0 "
                   " 0 0 1 ";

     for (int i=0; i<3; ++i)
       for (int j=0; j<3; ++j) 
         CPPUNIT_ASSERT_DOUBLES_EQUAL(Id[i][j], T[0][i][j], 1e-10);
         
   }

#undef __FUNC__
#define __FUNC__ "DistanceFromTargetTest::test_quad_with_default_target_matrix"
   void test_quad_with_default_target_matrix()
   {
     MsqError err;
     Matrix3D tri_m3d, quad_m3d, tet_m3d, hex_m3d; 
        
     TargetCalculator::initialize_default_target_matrices(tri_m3d, quad_m3d, tet_m3d, hex_m3d);

     quadPatch.allocate_corner_matrices(err); MSQ_CHKERR(err);
     MsqMeshEntity* elem = quadPatch.get_element_array(err);
     MsqTag* tag = elem[0].get_tag();
     tag->target_matrix(0) = quad_m3d;
     tag->target_matrix(1) = quad_m3d;
     tag->target_matrix(2) = quad_m3d;
     tag->target_matrix(3) = quad_m3d;

     Matrix3D T[MSQ_MAX_NUM_VERT_PER_ENT];
     double c_k[MSQ_MAX_NUM_VERT_PER_ENT];
     this->compute_T_matrices(elem[0], quadPatch, T, 4, c_k, err); MSQ_CHKERR(err);

//     cout << T[0] << endl;

     Matrix3D Id = " 1 0 0 "
                   " 0 1 0 "
                   " 0 0 1 ";

     for (int i=0; i<3; ++i)
       for (int j=0; j<3; ++j) 
         CPPUNIT_ASSERT_DOUBLES_EQUAL(Id[i][j], T[0][i][j], 1e-10);
         
   }
  
#undef __FUNC__
#define __FUNC__ "DistanceFromTargetTest::test_tet_with_default_target_matrix"
   void test_tet_with_default_target_matrix()
   {
     MsqError err;
     Matrix3D tri_m3d, quad_m3d, tet_m3d, hex_m3d; 
        
     TargetCalculator::initialize_default_target_matrices(tri_m3d, quad_m3d, tet_m3d, hex_m3d);

     tetPatch.allocate_corner_matrices(err); MSQ_CHKERR(err);
     MsqMeshEntity* elem = tetPatch.get_element_array(err);
     MsqTag* tag = elem[0].get_tag();
     tag->target_matrix(0) = tet_m3d;
     tag->target_matrix(1) = tet_m3d;
     tag->target_matrix(2) = tet_m3d;
     tag->target_matrix(3) = tet_m3d;

     double coords[3];
     MsqVertex* vert_array = tetPatch.get_vertex_array(err);
     double fac = pow(sqrt(2.), 1./3.);
     coords[0] = 0; coords[1] = 0; coords[2] = 0;
     vert_array[0] = coords;
     coords[0] = fac; coords[1] = 0; coords[2] = 0;
     vert_array[1] = coords;
     coords[0] = fac*.5; coords[1] = fac*sqrt(3.)/2.; coords[2] = 0;
     vert_array[2] = coords;
     coords[0] = fac*1./2.; coords[1] = fac*sqrt(3.)/6.; coords[2] = fac*sqrt(2.)/sqrt(3.);
     vert_array[3] = coords;

     Matrix3D T[MSQ_MAX_NUM_VERT_PER_ENT];
     double c_k[MSQ_MAX_NUM_VERT_PER_ENT];
     this->compute_T_matrices(elem[0], tetPatch, T, 4, c_k, err); MSQ_CHKERR(err);

//     cout << T[0] << endl;

     Matrix3D Id = " 1 0 0 "
                   " 0 1 0 "
                   " 0 0 1 ";

     for (int i=0; i<3; ++i)
       for (int j=0; j<3; ++j) 
         CPPUNIT_ASSERT_DOUBLES_EQUAL(Id[i][j], T[0][i][j], 1e-10);
         
   }

#undef __FUNC__
#define __FUNC__ "DistanceFromTargetTest::test_compute_T_matrices"
   void test_compute_T_matrices()
   {
     MsqError err;
     ShapeGuides811 target_calc;
     Matrix3D T[MSQ_MAX_NUM_VERT_PER_ENT];
     double c_k[MSQ_MAX_NUM_VERT_PER_ENT];

     
     
     // W is a singular matrix in 2D ... code needs fixing.
     MsqMeshEntity* elem = oneTri.get_element_array(err); MSQ_CHKERR(err);
     target_calc.compute_target_matrices(oneTri, err); MSQ_CHKERR(err);
     this->compute_T_matrices(elem[0], oneTri, T, 3, c_k, err); MSQ_CHKERR(err);
     
     Matrix3D T_check[3];

     T_check[0] = "  1.90637    0        -0.169031 "
                  "  0.953184  -2.75161  -0.507093 "
                  " -0.953184   1.65096  -0.845154 ";

     T_check[1] = " -0.953184 -1.65096  -0.169031 "
                  " -2.85955   0.550321 -0.507093 "
                  "  1.90637   0        -0.845154 ";

     T_check[2] = " -0.953184  1.65096  -0.169031 "
                  "  1.90637   2.20128  -0.507093 "
                  " -0.953184 -1.65096  -0.845154 ";

//       oneTri.print();
//       cout << T[0] << endl << T[1] << endl << T[2] << endl;
     
     for (int t=0; t<3; ++t) 
       for (int i=0; i<3; ++i)
         for (int j=0; j<3; ++j) {
//            cout << t << " " << i << " " << j << " " << T_check[t][i][j] << " " << T[t][i][j] << endl;
           CPPUNIT_ASSERT_DOUBLES_EQUAL(T_check[t][i][j], T[t][i][j], .00001);
         }

     double value1;     
     I_DFT mu1;
     mu1.evaluate_element(oneTri, &elem[0], value1, err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(3.17393, value1, .00001);
     
     double value2;     
     sI_DFT mu2;
     mu2.evaluate_element(oneTri, &elem[0], value2, err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(4.43943, value2, .00001);
     
     double value3;     
     RI_DFT mu3;
     mu3.evaluate_element(oneTri, &elem[0], value3, err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(6.20141, value3, .00001);
     
     double value4;
     sRI_DFT mu4;
     mu4.evaluate_element(oneTri, &elem[0], value4, err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(4.07917, value4, .00001);

   }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(DistanceFromTargetTest, "DistanceFromTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(DistanceFromTargetTest, "Unit");
