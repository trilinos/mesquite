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

#include "Mesquite.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "MeshImpl.hpp"

#include "cppunit/extensions/HelperMacros.h"

#include <iostream>

class ExodusTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(ExodusTest);
   CPPUNIT_TEST (test_elements);
   CPPUNIT_TEST_SUITE_END();

private:
  Mesquite::MeshImpl *mMesh;
  
public:
   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
#ifdef MSQ_USING_EXODUS
    Mesquite::MsqError err;
    
      // Read a Exodus Mesh file -- 10 triangles, 2 free vertices
    mMesh = new Mesquite::MeshImpl;
    mMesh->read_exodus("../../meshFiles/3D/CUBIT/test.g", err);
#endif
  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {
  }
  
public:
  ExodusTest()
    {}
  
  void test_elements()
  {
#ifdef MSQ_USING_EXODUS
    Mesquite::MsqError err;
    
      // Add mesh to a MeshSet.
    Mesquite::MeshSet mesh_set;
    mesh_set.add_mesh(mMesh, err); CPPUNIT_ASSERT(!err.errorOn);

      // Get the number of vertices
    std::cout << "Number of vertices: "
              << mMesh->get_total_vertex_count(err) << std::endl;
      // Get the number of elements
    std::cout << "Number of elements: "
              << mMesh->get_total_element_count(err) << std::endl;
#endif
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ExodusTest, "ExodusTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ExodusTest, "Unit");
