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

class VtkTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(VtkTest);
   CPPUNIT_TEST (test_elements);
   CPPUNIT_TEST_SUITE_END();

private:
  Mesquite::MeshImpl *mMesh;
  
public:
   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
    Mesquite::MsqError err;
    
      // Read a Vtk Mesh file -- 10 triangles, 2 free vertices
    mMesh = new Mesquite::MeshImpl;
    mMesh->read_vtk("../../meshFiles/2D/VTK/equil_tri2.vtk", err);
  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {
  }
  
public:
  VtkTest()
    {}
  
  void test_elements()
  {
    Mesquite::MsqError err;
    
      // Add mesh to a MeshSet.
    Mesquite::MeshSet mesh_set;
    mesh_set.add_mesh(mMesh, err); CPPUNIT_ASSERT(!err.errorOn);
    
      // Retrieve a patch
    Mesquite::PatchData pd;
    Mesquite::PatchDataParameters pd_params;
    pd_params.set_patch_type(Mesquite::PatchData::ELEMENTS_ON_VERTEX_PATCH,
                             err, 1, 0);
    pd_params.add_culling_method(Mesquite::PatchData::NO_BOUNDARY_VTX);
    
    mesh_set.get_next_patch(pd, pd_params, err); CPPUNIT_ASSERT(!err.errorOn);
//    mesh_set.write_gnuplot("toto", err); CPPUNIT_ASSERT(!err.errorOn);
    
    int free_vtx = pd.num_free_vertices(err); CPPUNIT_ASSERT(!err.errorOn);
//    std::cout << "nb of free vertices: " << free_vtx << std::endl;
    CPPUNIT_ASSERT( free_vtx == 1 );
    
    Mesquite::MsqMeshEntity* element_array =  pd.get_element_array(err); CPPUNIT_ASSERT(!err.errorOn);
    size_t num_elements = pd.num_elements();
    CPPUNIT_ASSERT( num_elements == 6 );
    
    Mesquite::MsqVertex* vtx_array = pd.get_vertex_array(err); CPPUNIT_ASSERT(!err.errorOn);
    size_t num_vertices = pd.num_vertices();
    CPPUNIT_ASSERT( num_vertices == 7 );
    
    CPPUNIT_ASSERT( tri_check_validity(element_array, num_elements, vtx_array, num_vertices) == 1 );
    
    mesh_set.get_next_patch(pd, pd_params, err); CPPUNIT_ASSERT(!err.errorOn);
    
    element_array =  pd.get_element_array(err); CPPUNIT_ASSERT(!err.errorOn);
    num_elements = pd.num_elements();
    CPPUNIT_ASSERT( num_elements == 6 );
    
    vtx_array = pd.get_vertex_array(err); CPPUNIT_ASSERT(!err.errorOn);
    num_vertices = pd.num_vertices();
    CPPUNIT_ASSERT( num_vertices == 7 );
    
    CPPUNIT_ASSERT( tri_check_validity(element_array, num_elements, vtx_array, num_vertices) == 1 );
  }
  
  int tri_check_validity(Mesquite::MsqMeshEntity* element_array,
                         size_t num_elements,
                         Mesquite::MsqVertex* vtx_array,
                         size_t num_vertices)
  {
       
      /* check that the simplicial mesh is still valid, 
         based on right handedness. Returns a 1 or a 0 */
    int valid = 1;
    double dEps = 1.e-13;
    
    double x1, x2, x3, y1, y2, y3;// z1, z2, z3;
    std::vector<size_t> vertex_indices;
    
    for (size_t i=0;i<num_elements;i++)
    {
      element_array[i].get_vertex_indices(vertex_indices);
      
      x1 = vtx_array[vertex_indices[0]][0];
      y1 = vtx_array[vertex_indices[0]][1];
      x2 = vtx_array[vertex_indices[1]][0];
      y2 = vtx_array[vertex_indices[1]][1];
      x3 = vtx_array[vertex_indices[2]][0];
      y3 = vtx_array[vertex_indices[2]][1];
      
      double a = x2*y3 - x3*y2;
      double b = y2 - y3;
      double c = x3 - x2;
      
      double area = .5*(a+b*x1+c*y1);
      if (area < dEps) {
          //          printf("x1 y1 = %f %f\n",x1,y1);
          //          printf("x2 y3 = %f %f\n",x2,y2);
          //          printf("x3 y3 = %f %f\n",x3,y3);
          //          printf("area = %f\n",area);
        valid=0;
      }
    }
    
    return(valid);
  }
  
  int tet_validity_check(Mesquite::MsqMeshEntity* element_array,
                         size_t num_elements,
                         Mesquite::MsqVertex *vtx_array)
  {
    int valid = 1;
    double dEps = 1.e-13;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
    std::vector<size_t> vertex_indices;
    
    for (size_t i=0;i<num_elements;i++)
    {
      element_array[i].get_vertex_indices(vertex_indices);
      
      x1=vtx_array[vertex_indices[0]][0];
      y1=vtx_array[vertex_indices[0]][1];
      z1=vtx_array[vertex_indices[0]][2];
      
      x2=vtx_array[vertex_indices[1]][0];
      y2=vtx_array[vertex_indices[1]][1];
      z2=vtx_array[vertex_indices[1]][2];
      
      x3=vtx_array[vertex_indices[2]][0];
      y3=vtx_array[vertex_indices[2]][1];
      z3=vtx_array[vertex_indices[2]][2];
      
      x4=vtx_array[vertex_indices[3]][0];
      y4=vtx_array[vertex_indices[3]][1];
      z4=vtx_array[vertex_indices[3]][2];
      
      double dDX2 = x2 - x1;
      double dDX3 = x3 - x1;
      double dDX4 = x4 - x1;
      
      double dDY2 = y2 - y1;
      double dDY3 = y3 - y1;
      double dDY4 = y4 - y1;
      
      double dDZ2 = z2 - z1;
      double dDZ3 = z3 - z1;
      double dDZ4 = z4 - z1;
      
        /* dDet is proportional to the cell volume */
      double dDet = dDX2*dDY3*dDZ4 + dDX3*dDY4*dDZ2 + dDX4*dDY2*dDZ3
        - dDZ2*dDY3*dDX4 - dDZ3*dDY4*dDX2 - dDZ4*dDY2*dDX3 ;
      
        /* Compute a length scale based on edge lengths. */
      double dScale = ( sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) +
                             (z1-z2)*(z1-z2)) +
                        sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) +
                             (z1-z3)*(z1-z3)) +
                        sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) +
                             (z1-z4)*(z1-z4)) +
                        sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) +
                             (z2-z3)*(z2-z3)) +
                        sqrt((x2-x4)*(x2-x4) + (y2-y4)*(y2-y4) +
                             (z2-z4)*(z2-z4)) +
                        sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) +
                             (z3-z4)*(z3-z4)) ) / 6.;
      
        /* Use the length scale to get a better idea if the tet is flat or
           just really small. */
      if (fabs(dScale) < dEps)
      {
        return(valid = 0);
      }
      else
      {
        dDet /= (dScale*dScale*dScale);
      }
      
      if (dDet > dEps)
      {
        valid = 1;
      }
      else if (dDet < -dEps)
      {
        valid = -1;
      }
      else
      {
        valid = 0;
      }
    }  // end for i=1,numElements
    
    return(valid);
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VtkTest, "VtkTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VtkTest, "Unit");
