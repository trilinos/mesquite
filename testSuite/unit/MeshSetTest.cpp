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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 22-Nov-02 at 08:05:56
//  LAST-MOD:  3-Dec-02 at 15:24:27 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MeshSetTest.cpp

Unit testing of various functions in the MeshSet class. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "MeshImpl.hpp"

#include "cppunit/extensions/HelperMacros.h"

#include <list>
#include <iterator>

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else 
#include <iostream>
using std::cout;
using std::endl;
#endif

using namespace Mesquite;

class MeshSetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(MeshSetTest);
  CPPUNIT_TEST (test_add_multiple_meshes);
  CPPUNIT_TEST_SUITE_END();
  
private:
  
  Mesquite::MeshImpl *tri8;
  Mesquite::MeshImpl *quad4;
  
public:
    /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
    MsqPrintError err(cout);
    
      // Read a vtk file -- square meshed by 8 triangles
    tri8 = new Mesquite::MeshImpl;
    tri8->read_vtk("../../meshFiles/2D/VTK/square_tri_2.vtk", err);
    CPPUNIT_ASSERT(!err);
    
      // Read a vtk file
      // -- square meshed by 4 quads, adjacent to the previous mesh
    quad4 = new Mesquite::MeshImpl;
    quad4->read_vtk("../../meshFiles/2D/VTK/four_more_quads.vtk", err);
    CPPUNIT_ASSERT(!err);
  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {
  }
  
public:
  MeshSetTest()
   {}

   /*this function test the MeshSet concept of concatenating several meshes.*/
   void test_add_multiple_meshes()
   {
     MsqPrintError err(cout);
     
       /* Adds 2 adjacent meshes to the MeshSet. */
     MeshSet mesh_set;
     mesh_set.add_mesh(tri8, err); 
     CPPUNIT_ASSERT(!err);
     mesh_set.add_mesh(quad4, err); 
     CPPUNIT_ASSERT(!err);
     
       /* Retrieves a global patch */
     PatchData pd;
     PatchDataParameters pd_params;
     pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err, 0, 0);
     CPPUNIT_ASSERT(!err);
     mesh_set.get_next_patch(pd, pd_params, err);
     CPPUNIT_ASSERT(!err);
     
     int num_vtx = pd.num_vertices();
     CPPUNIT_ASSERT( num_vtx == 18 );
   }
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshSetTest, "MeshSetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshSetTest, "Unit");
