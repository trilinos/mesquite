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
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 18-Oct-04 by J.Kraftcheck
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifndef MSQ_USE_OLD_IO_HEADERS
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#endif

#ifdef MSQ_USE_OLD_C_HEADERS
#include <cstdlib>
#else
#include <stdlib.h>
#endif


#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"

// algorithms
#include "MeshTransform.hpp"
#include "Matrix3D.hpp"
#include "Vector3D.hpp"

using namespace Mesquite;


int main(int argc, char* argv[])
{
  Mesquite::MsqPrintError err(cout);

  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
    //mesh->read_exodus("transformed_mesh.exo", err);
  mesh->read_vtk("../../meshFiles/2D/VTK/tfi_horse10x4-12.vtk", err);
  if (err) return 1;
  
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err);
  if (err) return 1;
    //create the matrix for affine transformation
  double array_entries[9];
  array_entries[0]=0; array_entries[1]=1; array_entries[2]=0;
  array_entries[3]=1; array_entries[4]=0; array_entries[5]=0;
  array_entries[6]=0; array_entries[7]=0; array_entries[8]=1;
    //create the translation vector
  Matrix3D my_mat(array_entries);
  Vector3D my_vec(0, 0 , 10);
  MeshTransform my_transform(my_mat, my_vec, err);
    //my_transform.add_culling_method(PatchData::NO_INTERIOR_VTX);
    //my_transform.add_culling_method(PatchData::NO_BOUNDARY_VTX);
  if (err) return 1;
    //mesh->write_exodus("original_mesh.exo", err);
  my_transform.loop_over_mesh(mesh_set1, err);
    //mesh->write_exodus("transformed_mesh.exo", err);
  mesh->write_vtk("transformed_mesh.vtk", err);
}
