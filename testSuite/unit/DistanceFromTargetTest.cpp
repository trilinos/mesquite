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
#include "sI_DFT.hpp"

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
  CPPUNIT_TEST (test_compute_T_matrices);
  CPPUNIT_TEST_SUITE_END();
  
private:
  
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
    
     /* Our tet patch is made of two tets.  tet_1 is a perfect
        equilateral (the ideal for most metrics).  tet_2 is an arbitrary
        tet.
     */
    create_qm_two_tet_patch(tetPatch,err);MSQ_CHKERR(err);
    
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
#define __FUNC__ "DistanceFromTargetTest::test_compute_T_matrices"
   void test_compute_T_matrices()
   {
     MsqError err;
     DefaultTargetCalculator target_calc;
     Matrix3D T[MSQ_MAX_NUM_VERT_PER_ENT];
     double c_k[MSQ_MAX_NUM_VERT_PER_ENT];

     // W is a singular matrix in 2D ... code needs fixing.
     MsqMeshEntity* tri_elems = idealTri.get_element_array(err); MSQ_CHKERR(err);
     target_calc.compute_target_matrices(idealTri, err); MSQ_CHKERR(err);
     this->compute_T_matrices(tri_elems[0], idealTri, T, 3, c_k, err); MSQ_CHKERR(err);

     MsqMeshEntity* elems = idealTet.get_element_array(err); MSQ_CHKERR(err);
     target_calc.compute_target_matrices(idealTet, err); MSQ_CHKERR(err);
     this->compute_T_matrices(elems[0], idealTet, T, 4, c_k, err); MSQ_CHKERR(err);

     for (int i=0; i<4; ++i)
       cout << "T["<<i<<"]:\n" << T[i] << endl;
   }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(DistanceFromTargetTest, "DistanceFromTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(DistanceFromTargetTest, "Unit");
